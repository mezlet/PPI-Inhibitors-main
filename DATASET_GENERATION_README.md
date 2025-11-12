# Researcher Dataset Generation Script

This directory contains a script to generate the researcher's dataset with exact specifications for the PPI Inhibitors project.

## Dataset Specifications

The generated dataset contains:

- **Total Examples**: 11,127
- **Positive Examples** (Inhibitors): 714
- **Negative Examples** (Non-inhibitors): 10,413
- **Unique Protein Complexes**: 22
- **Unique Inhibitors**: 606 (from positive examples)

## Usage Options

There are two ways to generate the researcher's dataset:

### Option 1: Standalone Script (`generate_researcher_dataset.py`)

**Purpose**: A standalone Python script that can be run from the command line.

**Requirements**:
- Python 3.6+
- No external libraries required (uses only standard library)

**Usage**:
```bash
python3 generate_researcher_dataset.py
```

### Option 2: Integrated in Complete Pipeline Notebook

**Purpose**: The dataset generation is now integrated into the complete pipeline notebook for seamless workflow.

**Location**: `Complete_PPI_Inhibitors_Pipeline_End_To_End.ipynb` → Section 7.5

**Usage**:
1. Open the notebook in Google Colab or Jupyter
2. Navigate to Section 7.5: "Generate Researcher Dataset"
3. Run the cells in that section
4. The dataset will be automatically generated and verified

**Features**:
- Automatic detection of existing dataset (won't regenerate if already exists)
- Integrated verification after generation
- Ready to use immediately in subsequent training steps
- All functions prefixed with `_nb` to avoid naming conflicts

### Input Files

The script requires the following files in the `Data/` directory:

1. **2p2iInhibitorsSMILES.txt** - Contains positive examples (known inhibitors)
   - Format: `complex_name inhibited_complex pdb_id ligand_id smiles label`
   - Total: 956 entries from 35 unique complexes

2. **BindersWithComplexname.csv** - Contains binders used as negative examples
   - Format: `Complexname,Binders SMILES`
   - Total: 11,789 entries

### Output File

- **Data/researcher_dataset.txt** - The generated dataset
  - Format: `complex_name target_complex compound_id label`
  - Label: `1.0` for positive (inhibitor), `0.0` for negative (non-inhibitor)

### Algorithm

The script uses a multi-step approach to generate the dataset:

#### Step 1: Complex Selection
- Ranks all 35 complexes by number of unique inhibitors
- Selects the top 22 complexes to maximize unique inhibitor coverage
- Ensures exactly 22 unique protein complexes in the final dataset

#### Step 2: Positive Example Generation
1. **First Pass**: Selects examples to ensure exactly 606 unique inhibitors
   - One example per unique inhibitor
2. **Second Pass**: Adds additional examples from the 606 selected inhibitors
   - Reaches the target of 714 total positive examples
   - Maintains the 606 unique inhibitor constraint

#### Step 3: Negative Example Generation
Generates 10,413 negative examples using two strategies:

1. **Binders Strategy** (~75% of negatives)
   - Uses compounds that bind to proteins but don't inhibit the complex
   - Sourced from BindersWithComplexname.csv
   - Limited to the 22 selected complexes

2. **Cross-Complex Strategy** (~25% of negatives)
   - Pairs inhibitors with complexes they don't inhibit
   - Creates biologically plausible negative examples
   - Prevents the model from learning trivial patterns

#### Step 4: Dataset Combination
- Combines positive and negative examples
- Shuffles the dataset randomly
- Writes to output file

### Features

- **Reproducible**: Uses fixed random seed (42) for consistent results
- **Validated**: Automatically verifies all specifications are met
- **Comprehensive**: Includes detailed logging of each step
- **Efficient**: Optimized algorithms prevent infinite loops

### Verification

The script automatically verifies that the generated dataset meets all specifications:

```
============================================================
DATASET VERIFICATION REPORT
============================================================
TOTAL EXAMPLES: 11127 (target: 11127)
POSITIVE EXAMPLES: 714 (target: 714)
NEGATIVE EXAMPLES: 10413 (target: 10413)
UNIQUE PROTEIN COMPLEXES: 22 (target: 22)
UNIQUE INHIBITORS (from positive examples): 606 (target: 606)
============================================================

✓ ALL CHECKS PASSED!
```

### Top 22 Selected Complexes

The script selects these complexes (ordered by number of inhibitors):

1. 3UVW_A_2_B (206 inhibitors)
2. 2B4J_A_2_B (71 inhibitors)
3. 2RNY_A_2_B (61 inhibitors)
4. 4QC3_A_2_B (52 inhibitors)
5. 4QC3_A_2_C (52 inhibitors)
6. 1YCR_A_2_B (51 inhibitors)
7. 4AJY_C_2_B (30 inhibitors)
8. 4AJY_C_2_H (30 inhibitors)
9. 4AJY_C_2_V (30 inhibitors)
10. 1NW9_A_2_B (29 inhibitors)
11. 4GQ6_A_2_B (23 inhibitors)
12. 2E3K_A_2_B (23 inhibitors)
13. 2E3K_A_2_C (23 inhibitors)
14. 2E3K_A_2_Q (23 inhibitors)
15. 4YY6_A_2_Z (22 inhibitors)
16. 1BXL_A_2_B (21 inhibitors)
17. 4ESG_A_2_B (16 inhibitors)
18. 4ESG_A_2_D (16 inhibitors)
19. 2FLU_X_2_P (12 inhibitors)
20. 1YCQ_A_2_B (11 inhibitors)
21. 3D9T_A_2_C (10 inhibitors)
22. 2XA0_A_2_B (9 inhibitors)

These 22 complexes provide 615 unique inhibitors, from which 606 are selected for the dataset.

## Example Output

```
3UVW_A_2_B 3UVW_A_2_B WSH 1.0
1YCR_A_2_B 1YCR_A_2_B 2227 0.0
2E3K_A_2_B 2E3K_A_2_B X8U 0.0
2B4J_A_2_B 2B4J_A_2_B 3I1 1.0
...
```

## Notes

- The random seed is set to 42 for reproducibility
- The script handles edge cases like insufficient unique inhibitors
- Cross-complex negatives prevent the model from memorizing complex-compound pairs
- The 1:14.6 positive-to-negative ratio reflects the natural imbalance in drug discovery

## Citation

If you use this dataset generation script, please cite the original PPI Inhibitors paper:

```
[Citation information from research_paper.pdf]
```

## Author

Generated using Claude Code
Date: 2025-11-12

## License

[Same as parent repository]
