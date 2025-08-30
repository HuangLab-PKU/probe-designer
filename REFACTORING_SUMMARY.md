# Refactoring Summary: Thermal Filter and Sequence Handling

## Overview
This refactoring moves thermodynamic filtering functionality from `search_strategies.py` to `filtering.py` and ensures proper reverse complement sequence handling for probe design, with support for padlock probe arm sequences and enhanced specificity filtering.

## Key Changes

### 1. Thermal Filter Function (`src/filtering.py`)
- **New function**: `thermal_filter()` - Comprehensive thermodynamic screening function
- **Parameters**:
  - `arm_3prime`: 3' arm sequence
  - `arm_5prime`: 5' arm sequence
  - `sequence_type`: "DNA" or "RNA" (default: "DNA")
  - `target_type`: "DNA" or "RNA" (default: "RNA")
  - `target_sequence`: Target sequence for RNA secondary structure analysis
  - `**kwargs`: Optional parameter overrides
- **Returns**: Dictionary with filter results and thermodynamic parameters

### 2. Moved Functions from `search_strategies.py` to `filtering.py`
The following functions were moved and consolidated into `thermal_filter()`:
- `_pre_filter_sequence()`
- `_calculate_tm()`
- `_check_tm_constraints()`
- `_is_valid_binding_site_comprehensive()`
- `_calculate_tm_3_prime()`
- `_calculate_tm_5_prime()`
- `_check_rna_structure()`

### 3. Configuration Changes
- **Moved parameters** from `SearchConfig` to `FilterConfig`:
  - `min_g_content`, `max_g_content`
  - `max_consecutive_g`
  - `min_tm`, `max_tm`
  - `max_tm_diff`
  - `min_free_energy`
  - `check_rna_structure`

### 4. Padlock Probe Arm Sequence Handling
- **Arm sequence processing**: All strategies now generate separate 3' and 5' arm sequences
- **Correct orientation**: 3' arm is front half, 5' arm is back half (for padlock probes)
- **Dual storage**: Both arm sequences and full probe sequence are stored
- **Thermal analysis**: Separate Tm calculations for each arm

### 5. Enhanced Specificity Filtering
- **Renamed function**: `post_blast_filter()` → `specificity_filter()`
- **Dual criteria**: Checks both organism specificity and sequence complementarity
- **Frame-based complementarity**: Uses BLAST frame information (`frame[1] = -1`) to identify complementary strands
- **Enhanced validation**: More comprehensive specificity validation

### 6. Updated Search Strategies
All search strategies now:
- Extract target sequence from genome/transcript
- Generate probe sequence as reverse complement
- Split probe into 3' and 5' arms
- Apply thermal filter to arm sequences
- Store both arm sequences and full sequences in results

### 7. Configuration Files Updated
All configuration files (`config_*.json`) have been updated to:
- Remove thermodynamic parameters from `search` section
- Add thermodynamic parameters to `filter` section
- Maintain consistent structure across all config files

## Benefits

1. **Centralized filtering**: All thermodynamic filtering logic is now in one place
2. **Padlock probe support**: Proper handling of 3' and 5' arm sequences
3. **Enhanced specificity**: Dual criteria for organism and complementarity validation
4. **Frame-based detection**: Accurate complementarity detection using BLAST frame information
5. **Better separation of concerns**: Search strategies focus on finding positions, filtering handles thermodynamics
6. **Improved maintainability**: Changes to filtering logic only need to be made in one place
7. **English documentation**: All comments and documentation are now in English

## Usage Example

```python
from src.filtering import SequenceFilter
from src.config import FilterConfig, BlastConfig

# Create filter
filter_config = FilterConfig(
    min_g_content=0.3,
    max_g_content=0.7,
    min_tm=50.0,
    max_tm=80.0
)
blast_config = BlastConfig()
filter = SequenceFilter(filter_config, blast_config)

# Filter with arm sequences
arm_3prime = "ATCGATCGATCGATCGATCG"
arm_5prime = "GCTAGCTAGCTAGCTAGCTA"
target_seq = "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"

result = filter.thermal_filter(
    arm_3prime,
    arm_5prime,
    sequence_type="DNA",
    target_type="RNA",
    target_sequence=target_seq
)

if result['passed']:
    print(f"Sequence passed with Tm: {result['tm']:.1f}°C")
    print(f"3' arm Tm: {result['tm_3prime']:.1f}°C")
    print(f"5' arm Tm: {result['tm_5prime']:.1f}°C")
else:
    print(f"Failed checks: {result['failed_checks']}")
```

## Testing
- `test_specificity_filter.py`: Tests the specificity filter function and arm sequence handling
- `test_thermal_filter.py`: Tests the thermal filter function
- `test_reverse_complement.py`: Tests reverse complement sequence handling

## Key Features

### Padlock Probe Design
- **3' arm**: Front half of probe sequence
- **5' arm**: Back half of probe sequence
- **Thermal analysis**: Separate Tm calculations for each arm
- **RNA structure**: Analysis based on target sequence

### Specificity Validation
- **Organism specificity**: Checks against target organisms
- **Complementarity**: Verifies probe-reference complementarity using BLAST frame information
- **Frame-based detection**: Uses `frame[1] = -1` to identify complementary strands
- **Gene specificity**: Checks for non-target gene alignments within target organisms
- **Enhanced BLAST integration**: Improved alignment analysis with frame parsing

### BLAST Species Configuration
- **Configurable species**: BLAST search target species can be set in config or command line
- **Flexible targeting**: Support single species, multiple species, or automatic species detection
- **Priority system**: BLAST species > filter target_organisms > config species
- **Command line override**: `--blast_species` parameter for runtime configuration
- **No alignment count limits**: Removed max_alignments restriction for better flexibility
