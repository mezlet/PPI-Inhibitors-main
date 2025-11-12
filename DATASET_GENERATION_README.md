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

The script requires the following file in the `Data/` directory:

**WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt** - Precomputed dataset containing all examples
   - Format: `complex_name target_complex compound_id label`
   - Note: Compound IDs can contain spaces (e.g., "quingestanol acetate")
   - Total: 15,695 examples
     - 857 positives (with 606 unique inhibitors)
     - 14,838 negatives
     - 22 unique protein complexes

### Output File

- **Data/researcher_dataset.txt** - The generated dataset
  - Format: `complex_name target_complex compound_id label`
  - Label: `1.0` for positive (inhibitor), `0.0` for negative (non-inhibitor)

### Algorithm

The script uses an intelligent sampling approach from the precomputed dataset:

#### Step 1: Read Precomputed Dataset
- Loads `WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt`
- Separates into positives (857) and negatives (14,838)
- Verifies 22 unique complexes and 606 unique inhibitors are present
- Handles compound IDs with spaces correctly

#### Step 2: Sample Positive Examples
Samples 714 positive examples from 857 available while maintaining all 606 unique inhibitors:

1. **First Pass**: Ensures all unique inhibitors are represented
   - Groups examples by inhibitor
   - Randomly selects one example for each of the 606 unique inhibitors

2. **Second Pass**: Adds remaining examples to reach 714 total
   - Samples from remaining examples
   - Can include multiple examples for some inhibitors
   - Maintains diversity across complexes

#### Step 3: Sample Negative Examples
- Randomly samples 10,413 negatives from 14,838 available
- Maintains representation across all 22 complexes
- Preserves the mix of binders and cross-complex negatives from precomputed data

#### Step 4: Dataset Combination
- Combines 714 positives + 10,413 negatives = 11,127 total
- Shuffles the dataset randomly (using seed 42 for reproducibility)
- Writes to output file

### Features

- **Uses Precomputed Data**: Samples from the actual research dataset
- **Reproducible**: Uses fixed random seed (42) for consistent results
- **Validated**: Automatically verifies all specifications are met
- **Comprehensive**: Includes detailed logging of each step
- **Fast**: Direct sampling is much faster than regenerating from raw sources
- **Accurate**: Uses the exact same data that was precomputed for research

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

### The 22 Protein Complexes

The precomputed dataset contains these 22 protein complexes (ordered by unique inhibitor count):

1. 3UVW_A_2_B (201 unique inhibitors, 201 examples)
2. 2B4J_A_2_B (65 unique inhibitors, 65 examples)
3. 2RNY_A_2_B (61 unique inhibitors, 61 examples)
4. 4QC3_A_2_C (52 unique inhibitors, 104 examples)
5. 1YCR_A_2_B (51 unique inhibitors, 51 examples)
6. 4AJY_C_2_B (30 unique inhibitors, 90 examples)
7. 4GQ6_A_2_B (23 unique inhibitors, 23 examples)
8. 2E3K_A_2_Q (22 unique inhibitors, 66 examples)
9. 4YY6_A_2_Z (22 unique inhibitors, 22 examples)
10. 4ESG_A_2_D (15 unique inhibitors, 30 examples)
11. 1NW9_A_2_B (13 unique inhibitors, 13 examples)
12. 2FLU_X_2_P (12 unique inhibitors, 12 examples)
13. 1BXL_A_2_B (12 unique inhibitors, 12 examples)
14. 1YCQ_A_2_B (11 unique inhibitors, 11 examples)
15. 3D9T_A_2_D (10 unique inhibitors, 28 examples)
16. 3WN7_A_2_M (9 unique inhibitors, 27 examples)
17. 1Z92_A_2_B (7 unique inhibitors, 7 examples)
18. 3TDU_A_2_F (5 unique inhibitors, 20 examples)
19. 3DAB_A_2_B (5 unique inhibitors, 5 examples)
20. 1F47_A_2_B (4 unique inhibitors, 4 examples)
21. 2XA0_A_2_B (3 unique inhibitors, 3 examples)
22. 1BKD_S_2_R (2 unique inhibitors, 2 examples)

These 22 complexes collectively contain 606 unique inhibitors across 857 positive examples in the precomputed dataset.

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
