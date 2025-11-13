# PPI Inhibitors Dataset Generation

This directory contains code to generate the training dataset used in the research paper:
**"Predicting small-molecule inhibition of protein complexes"**

## Overview

The dataset generation process creates a comprehensive dataset matching the specifications in the research paper, consisting of:

- **714 positive examples** from the 2P2I database (22 protein complexes, 608 unique inhibitors)
- **~10,413 negative examples** generated using three strategies
- **Feature extraction** including ligand fingerprints, protein sequence features, and interface features

## Dataset Composition

### Positive Examples (714 total)
- Source: 2P2I v2 database
- Experimentally verified inhibitors of protein complexes
- Filtered to remove complexes with only one inhibitor
- Each example: protein complex + compound known to inhibit it

### Negative Examples (~10,413 total)

#### Strategy 1: Random Pairing (~857 examples)
- Random pairing of 2P2I complexes with compounds from 2P2I and SuperDRUG2
- Ensures selected compound is not a known inhibitor of that complex

#### Strategy 2: DBD5 Pairing (~1,714 examples)
- 2P2I compounds paired with DBD5 benchmark database complexes
- Uses 282 complexes from DBD5 (version 5.5)

#### Strategy 3: Binders that are NOT Inhibitors (~7,842 examples)
- Compounds that bind to proteins but don't inhibit the complex
- Sourced from Binding-DB database
- BLASTp search at >90% sequence identity
- Binding affinity: Ki, Kd, IC50 < 7.6 nM
- Tanimoto coefficient < 0.85 with known inhibitors

## Feature Extraction

The script extracts multiple types of features for each example:

### 1. Ligand Features (2,048 dimensions)
- **ECFP (Extended-Connectivity Fingerprint)** / Morgan Fingerprints
- Radius: 2 bonds
- Captures structural information of compounds

### 2. Protein Sequence Features (69 dimensions)
- **Amino Acid Composition (AAC)**: 20 dimensions
  - Frequency of each amino acid
- **Grouped k-mer composition**: 49 dimensions (k=2)
  - Groups amino acids by physicochemical properties (7 groups)
  - Features averaged across all chains in complex

### 3. Interface Features (211 dimensions)
- Frequency of amino acid pairs at protein-protein interface
- Interface defined as residues within 8Å distance
- 21×21 matrix (20 standard amino acids + unknown)

### Total Feature Dimension
- **2,328 dimensions** (pre-GNN features)
- **Note**: GNN features (512 dim) are computed during model training from 3D structure

## Prerequisites

```bash
# Install required packages
pip install rdkit-pypi
pip install biopython
pip install pandas numpy openpyxl
```

## Usage

### Basic Usage (Generate dataset only)

```bash
python generate_training_dataset.py
```

This will:
1. Load positive examples from 2P2I
2. Generate negative examples using all three strategies
3. Save dataset files to `./generated_data/`

### Advanced Usage

```bash
# Specify custom directories
python generate_training_dataset.py \
    --data_dir ./Data \
    --output_dir ./my_dataset \
    --seed 42

# Generate dataset without computing features (faster)
python generate_training_dataset.py --no_features

# Full feature extraction (slower, but includes all features)
python generate_training_dataset.py --data_dir ./Data --output_dir ./generated_data
```

### Command-line Arguments

- `--data_dir`: Directory containing input data files (default: `./Data`)
- `--output_dir`: Directory to save generated dataset (default: `./generated_data`)
- `--no_features`: Skip feature extraction, only generate dataset files
- `--seed`: Random seed for reproducibility (default: 42)

## Output Files

The script generates the following files in the output directory:

### Dataset Files
- **`dataset_all_examples.txt`**: Paper format (ComplexName TargetComplexName InhibitorName Label)
- **`dataset_complete.csv`**: Detailed CSV with all metadata
- **`dataset_positive.csv`**: Only positive examples
- **`dataset_negative.csv`**: Only negative examples

### Feature Files (if `--no_features` not specified)
- **`features.npy`**: Numpy array of shape (N, 2328) containing extracted features
- **`labels.npy`**: Numpy array of shape (N,) containing labels (0/1)

## Data Format

### Text File Format (dataset_all_examples.txt)
```
ComplexName TargetComplexName InhibitorCode Label
3UVW_A_2_B 3UVW_A_2_B WSH 1
1YCR_A_2_B 1YCR_A_2_B RANDOM_1 0
```

### CSV Format (dataset_complete.csv)
Columns:
- `complex_pair`: Full complex identifier
- `complex_name`: Base complex name
- `inhibitor_name`: Inhibitor identifier
- `inhibitor_code`: Short inhibitor code
- `smiles`: SMILES representation of compound
- `label`: 0 (negative) or 1 (positive)
- `strategy`: Generation strategy (for negative examples)

## Required Input Data Files

The script expects the following files in the data directory:

### Essential Files
- `2p2iInhibitorsSMILES.txt`: Positive examples with SMILES
- `2p2iComplexPairs.txt`: Protein complex sequences
- `BindersWithComplexname.csv` or `Binders With Tanimoto Similarity 0.85.csv`: Binder compounds

### Optional Files
- `approved_drugs_chemical_structure_identifiers.xlsx`: SuperDRUG2 compounds
- `DBD5_seq_dict`: DBD5 complex information
- `Pdb/*.pdb`: PDB structures for interface feature extraction

## Example: Loading Generated Dataset

```python
import pandas as pd
import numpy as np

# Load dataset
df = pd.read_csv('generated_data/dataset_complete.csv')
print(f"Total examples: {len(df)}")
print(f"Positive: {len(df[df['label']==1])}")
print(f"Negative: {len(df[df['label']==0])}")

# Load features
features = np.load('generated_data/features.npy')
labels = np.load('generated_data/labels.npy')

print(f"Features shape: {features.shape}")
print(f"Labels shape: {labels.shape}")
```

## Dataset Statistics (Expected)

```
Positive examples: ~714
Negative examples: ~10,413
  - Strategy 1 (Random pairing): ~857
  - Strategy 2 (DBD5 pairing): ~1,714
  - Strategy 3 (Binders): ~7,842
Total examples: ~11,127
Class balance: 714 / 10,413 (~1:14.6)
```

## Cross-Validation Strategy

The paper uses **Leave-One-Complex-Out (LOCO)** cross-validation:
- All examples of a complex are held out for testing
- Model trained on all other complexes
- Repeated for each of the 22 complexes

## Notes

1. **Class Imbalance**: The dataset has ~14x more negative than positive examples, reflecting real-world scenarios. The paper addresses this with weighted training.

2. **Hard Negative Examples**: Strategy 3 provides "hard" negatives (binders that don't inhibit), making the task more realistic and challenging.

3. **Feature Computation Time**:
   - Dataset generation only: ~1-2 minutes
   - With feature extraction: ~10-30 minutes (depending on PDB availability)

4. **PDB Structures**: Interface features require 3D structures. If PDB files are not available, interface features will be zero vectors.

5. **GNN Features**: The 512-dimensional GNN features are computed during model training using the graph neural network, not in this preprocessing step.

## Troubleshooting

### Issue: Missing data files
**Solution**: Download required files from the GitHub repository or paper supplementary materials.

### Issue: RDKit import error
**Solution**: Install RDKit: `pip install rdkit-pypi`

### Issue: Very few Strategy 3 negatives
**Solution**: Ensure `BindersWithComplexname.csv` or `Binders With Tanimoto Similarity 0.85.csv` is present in the Data directory.

### Issue: All interface features are zeros
**Solution**: Ensure PDB structures are available in `Data/Pdb/` directory. If not available, this is expected and won't affect other features.

## Citation

If you use this dataset generation code, please cite the paper:

```bibtex
@article{yaseen2024predicting,
  title={Predicting small-molecule inhibition of protein complexes},
  author={Yaseen, Adiba and Roy, Soumyadip and Akhter, Naeem and Ben-Hur, Asa and Minhas, Fayyaz},
  journal={Bioinformatics},
  year={2024}
}
```

## Contact

For questions or issues, please open an issue on the GitHub repository:
https://github.com/adibayaseen/PPI-Inhibitors

## License

This code is provided for research purposes. Please refer to the main repository for license information.
