#!/usr/bin/env python3
"""
DNA探针设计软件测试脚本
用于验证重构后的代码功能
"""

import os
import sys
import tempfile
import pandas as pd
from pathlib import Path

# 添加src目录到Python路径
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from config import ConfigManager
from database import DatabaseInterface
from search_strategies import BindingSiteSearcher
from filtering import SequenceFilter
from probe_assembly import ProbeAssembler, ProbeValidator
from utils import (
    create_output_directory, load_gene_list, save_gene_list,
    log_message, print_dependency_status, create_config_template,
    calculate_statistics, save_statistics
)


def test_config_manager():
    """测试配置管理器"""
    print("测试配置管理器...")
    
    # 创建默认配置
    config = ConfigManager()
    assert config.database.organism == "mouse"
    assert config.search.binding_site_length == 40
    assert config.filter.min_g_content == 0.4
    
    # 测试配置验证
    errors = config.validate_config()
    assert len(errors) == 0, f"配置验证失败: {errors}"
    
    print("✓ 配置管理器测试通过")


def test_utils():
    """测试工具函数"""
    print("测试工具函数...")
    
    # 测试基因列表加载和保存
    test_genes = ["CD3D", "CD4", "CD8A"]
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        for gene in test_genes:
            f.write(f"{gene}\n")
        temp_file = f.name
    
    try:
        loaded_genes = load_gene_list(temp_file)
        assert loaded_genes == test_genes
        print("✓ 基因列表加载测试通过")
    finally:
        os.unlink(temp_file)
    
    # 测试输出目录创建
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = create_output_directory(temp_dir, create_timestamp=False)
        assert os.path.exists(output_dir)
        print("✓ 输出目录创建测试通过")


def test_search_strategies():
    """测试搜索策略"""
    print("测试搜索策略...")
    
    from config import SearchConfig, FilterConfig
    
    search_config = SearchConfig(
        binding_site_length=20,
        max_binding_sites=5,
        search_strategy="brute_force"
    )
    
    filter_config = FilterConfig(
        min_g_content=0.3,
        max_g_content=0.7,
        min_tm=40.0,
        max_tm=70.0
    )
    
    # 创建测试序列
    test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    
    # 测试暴力搜索策略
    searcher = BindingSiteSearcher(search_config, filter_config)
    strategy = searcher.create_strategy("brute_force")
    
    binding_sites = strategy.search_binding_sites(test_sequence, "TEST_GENE")
    
    assert len(binding_sites) <= search_config.max_binding_sites
    if binding_sites:
        site = binding_sites[0]
        assert len(site['sequence']) == search_config.binding_site_length
        assert 'g_content' in site
        assert 'tm' in site
    
    print("✓ 搜索策略测试通过")


def test_probe_assembly():
    """测试探针组装"""
    print("测试探针组装...")
    
    from config import ProbeConfig
    
    probe_config = ProbeConfig(panel_type="PRISM")
    
    # 创建测试结合位点
    test_binding_sites = {
        "GENE1": [{
            'sequence': 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG',
            'position': 100,
            'g_content': 0.5,
            'tm': 55.0
        }],
        "GENE2": [{
            'sequence': 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA',
            'position': 200,
            'g_content': 0.6,
            'tm': 60.0
        }]
    }
    
    # 测试探针组装
    assembler = ProbeAssembler(probe_config)
    probe_df = assembler.assemble_probes(test_binding_sites)
    
    assert len(probe_df) == 2
    assert 'probe' in probe_df.columns
    assert 'gene' in probe_df.columns
    
    print("✓ 探针组装测试通过")


def test_probe_validation():
    """测试探针验证"""
    print("测试探针验证...")
    
    validator = ProbeValidator()
    
    # 测试有效序列
    valid_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    result = validator.validate_probe_sequence(valid_sequence)
    assert result['is_valid']
    
    # 测试无效序列（包含无效字符）
    invalid_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGX"
    result = validator.validate_probe_sequence(invalid_sequence)
    assert not result['is_valid']
    
    print("✓ 探针验证测试通过")


def test_config_template():
    """测试配置模板创建"""
    print("测试配置模板创建...")
    
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        create_config_template(temp_file)
        assert os.path.exists(temp_file)
        
        # 验证JSON格式
        config = ConfigManager(temp_file)
        assert config.database.organism == "mouse"
        assert config.search.binding_site_length == 40
        
        print("✓ 配置模板创建测试通过")
    finally:
        os.unlink(temp_file)


def test_dependency_check():
    """测试依赖检查"""
    print("测试依赖检查...")
    
    deps = check_dependencies()
    assert isinstance(deps, dict)
    assert 'pandas' in deps
    assert 'numpy' in deps
    
    print("✓ 依赖检查测试通过")


def create_test_gene_list():
    """创建测试基因列表"""
    test_genes = [
        "CD3D", "CD4", "CD8A", "CD19", "CD20",
        "CD14", "CD16", "CD56", "CD11c", "CD123"
    ]
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        for gene in test_genes:
            f.write(f"{gene}\n")
        return f.name


def run_integration_test():
    """运行集成测试"""
    print("运行集成测试...")
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as temp_dir:
        # 创建测试基因列表
        gene_list_file = create_test_gene_list()
        
        try:
            # 创建配置
            config = ConfigManager()
            config.database.organism = "mouse"
            config.search.max_binding_sites = 3  # 减少数量以加快测试
            config.output.output_dir = temp_dir
            config.output.save_intermediate = False  # 不保存中间文件
            
            # 创建输出目录
            output_dir = create_output_directory(temp_dir, create_timestamp=False)
            
            # 测试基因列表加载
            gene_list = load_gene_list(gene_list_file)
            assert len(gene_list) == 10
            
            # 测试配置验证
            errors = config.validate_config()
            assert len(errors) == 0
            
            print("✓ 集成测试通过")
            
        finally:
            os.unlink(gene_list_file)


def main():
    """主测试函数"""
    print("开始DNA探针设计软件测试...\n")
    
    try:
        # 运行各个测试
        test_config_manager()
        test_utils()
        test_search_strategies()
        test_probe_assembly()
        test_probe_validation()
        test_config_template()
        test_dependency_check()
        run_integration_test()
        
        print("\n🎉 所有测试通过！")
        print("重构后的代码功能正常。")
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
