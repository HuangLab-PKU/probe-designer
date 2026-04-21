"""TDD for probe_designer.config.loader.load_yaml_with_env.

Phase 1 will add this loader that expands ${VAR} and ${VAR:-default}
placeholders in YAML before parsing. Fail-fast on unresolved vars.

Expected signature:
    load_yaml_with_env(path: str | Path, env: Mapping[str, str] | None = None) -> dict

Semantics:
  - Reads YAML file as text
  - Expands ${VAR} using env (defaults to os.environ)
  - Expands ${VAR:-default} using env with fallback literal
  - Raises ConfigError with a pointer to .env.example when a ${VAR} is
    unresolved and has no default
  - Only strings are substituted; non-string YAML values untouched
"""
from __future__ import annotations

import pytest

from probe_designer.config.loader import load_yaml_with_env
from probe_designer.utils.errors import ConfigError


class TestBasicYamlLoading:
    def test_loads_plain_yaml_without_expansion(self, tmp_path):
        yaml_path = tmp_path / "plain.yaml"
        yaml_path.write_text("foo: 1\nbar: hello\n")
        result = load_yaml_with_env(yaml_path)
        assert result == {"foo": 1, "bar": "hello"}

    def test_returns_dict_for_nested_structure(self, tmp_path):
        yaml_path = tmp_path / "nested.yaml"
        yaml_path.write_text(
            "database:\n  host: localhost\n  port: 5432\n"
        )
        result = load_yaml_with_env(yaml_path)
        assert result["database"]["host"] == "localhost"
        assert result["database"]["port"] == 5432


class TestEnvVarExpansion:
    def test_plain_env_var_expanded(self, tmp_path):
        yaml_path = tmp_path / "env.yaml"
        yaml_path.write_text("email: ${USER_EMAIL}\n")
        result = load_yaml_with_env(yaml_path, env={"USER_EMAIL": "a@b.com"})
        assert result["email"] == "a@b.com"

    def test_default_used_when_var_missing(self, tmp_path):
        yaml_path = tmp_path / "default.yaml"
        yaml_path.write_text("level: ${LOG_LEVEL:-INFO}\n")
        result = load_yaml_with_env(yaml_path, env={})
        assert result["level"] == "INFO"

    def test_default_overridden_by_env(self, tmp_path):
        yaml_path = tmp_path / "override.yaml"
        yaml_path.write_text("level: ${LOG_LEVEL:-INFO}\n")
        result = load_yaml_with_env(yaml_path, env={"LOG_LEVEL": "DEBUG"})
        assert result["level"] == "DEBUG"

    def test_multiple_expansions_per_file(self, tmp_path):
        yaml_path = tmp_path / "multi.yaml"
        yaml_path.write_text(
            "database:\n  email: ${EMAIL}\n  key: ${KEY}\n  port: ${PORT:-5432}\n"
        )
        result = load_yaml_with_env(
            yaml_path, env={"EMAIL": "x@y.com", "KEY": "abc123"}
        )
        assert result["database"]["email"] == "x@y.com"
        assert result["database"]["key"] == "abc123"
        # Expanded "5432" text, post-YAML-parse, becomes int via YAML type inference
        assert result["database"]["port"] == 5432


class TestFailFast:
    def test_missing_var_no_default_raises_config_error(self, tmp_path):
        yaml_path = tmp_path / "missing.yaml"
        yaml_path.write_text("email: ${NONEXISTENT_VAR}\n")
        with pytest.raises(ConfigError):
            load_yaml_with_env(yaml_path, env={})

    def test_error_message_points_to_env_example(self, tmp_path):
        yaml_path = tmp_path / "missing.yaml"
        yaml_path.write_text("email: ${NONEXISTENT_VAR}\n")
        with pytest.raises(ConfigError, match=r"\.env\.example|environment"):
            load_yaml_with_env(yaml_path, env={})


class TestNonStringValuesUntouched:
    def test_numeric_values_preserved(self, tmp_path):
        yaml_path = tmp_path / "types.yaml"
        yaml_path.write_text(
            "batch_size: 100\n"
            "evalue: 1.0e-5\n"  # explicit decimal for PyYAML 1.1 float resolver
            "enabled: true\n"
            "items:\n  - a\n  - b\n"
        )
        result = load_yaml_with_env(yaml_path, env={})
        assert result["batch_size"] == 100
        assert result["evalue"] == 1.0e-5
        assert result["enabled"] is True
        assert result["items"] == ["a", "b"]


class TestPartialExpansionInString:
    def test_var_embedded_in_string_expanded(self, tmp_path):
        yaml_path = tmp_path / "embed.yaml"
        yaml_path.write_text("path: /data/${USER}/output\n")
        result = load_yaml_with_env(yaml_path, env={"USER": "alice"})
        assert result["path"] == "/data/alice/output"
