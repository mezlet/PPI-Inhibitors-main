# PPI Inhibitors Dataset Preprocessing Guide

## Overview

This guide explains how to preprocess the PPI Inhibitors dataset to match the exact specifications described in the research paper:

**"Predicting small-molecule inhibition of protein complexes"**
*Yaseen et al., 2024*

## Dataset Composition (From Research Paper)

### Training Dataset

The paper uses **714 positive examples (inhibitors)** from **22 protein complexes** from the 2P2I database.

| Component | Count | Description |
|-----------|-------|-------------|
| **Protein Complexes** | 22 | Protein-protein interaction complexes |
| **Positive Examples** | 714 | Experimentally verified inhibitors |
| **Negative Examples** | ~14,838 | Generated using 3 strategies |
| **Total Training** | ~15,552 | Used for Leave-One-Complex-Out validation |

### 22 Protein Complexes Used

From Table 2 of the paper:

1. **3DAB** - MDM4/P53
2. **3WN7** - MKEAP1/MNRF2
3. **2FLU** - KEAP1/NRF2
4. **1BKD** - HRAS/SOS1
5. **1YCQ** - XDM2/P53
6. **4ESG** - WDR5/MLL1
7. **4QC3** - BAZ2B/H4
8. **3TDU** - DCN1/UBC12
9. **1F47** - ZIPA/FTSZ
10. **2E3K** - BRD2-2/H4
11. **4AJY** - VHL/HIF1A
12. **3D9T** - CIAP1-BIR3/CASPASE-9
13. **2RNY** - CREBBP/H4
14. **3UVW** - BRD4-1/H4 (largest: 201 inhibitors)
15. **4YY6** - BRD9/H4
16. **1YCR** - MDM2/P53
17. **1BXL** - BCLXL/BAK
18. **2B4J** - INTEGRASE/LEDGF
19. **2XA0** - BCL2/BAX
20. **1Z92** - IL-2/IL-2R
21. **1NW9** - XIAP-BIR3/SMAC
22. **4GQ6** - MENIN/MLL

### Negative Example Generation (3 Strategies)

The paper uses **three strategies** to generate hard negative examples:

#### Strategy 1: Random Pairing (2P2I + SuperDRUG2)
- **Count**: ~857 examples
- **Method**: Randomly pair complexes from 2P2I with compounds from 2P2I and SuperDRUG2
- **Constraint**: Selected compound is NOT a known inhibitor of that complex

#### Strategy 2: DBD5 Complexes
- **Count**: ~1,714 examples
- **Method**: Randomly pair compounds from 2P2I with complexes from DBD5 benchmark database (v5.5)
- **Purpose**: Complexes where the compound is unlikely to be a binder

#### Strategy 3: Binders That Are Not Inhibitors
- **Count**: ~11,789 examples (after filtering)
- **Method**: Use BindingDB to find compounds that bind protein chains but are NOT inhibitors
- **Filtering Steps**:
  1. BLASTp search at >90% sequence identity
  2. Select binders with Ki, Kd, or IC50 < 7.6 nM
  3. Exclude compounds with Tanimoto coefficient ≥ 0.85 with known inhibitors
- **Purpose**: Create "hard" negatives that can bind but don't inhibit

### External Test Sets

#### External Test Set 1: Recent Publications
- **Count**: 28 inhibitors
- **Source**: Collected from recent literature
- **Purpose**: Test on novel proteins with low sequence similarity to training data

#### External Test Set 2: SARS-CoV-2 Spike/ACE2
- **Count**: 25 inhibitors
- **Complex**: SARS-CoV-2 Spike protein and Human ACE2
- **Source**: Hanson et al., 2020
- **Purpose**: Test on COVID-19 related protein-protein interaction

---

## Current Dataset Analysis

### Main Dataset File

**File**: `Data/WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt`

**Format**:
```
<Complex_Name> <Target_Complex> <Inhibitor_SMILES> <Label>
```

**Example**:
```
3UVW_A_2_B 3UVW_A_2_B WSH 1.0
3UVW_A_2_B 3UVW_A_2_B c1ccc(Br)cc1 0.0
```

### Analysis Results

After running the preprocessing script, we found:

| Metric | Paper | Actual File | Difference |
|--------|-------|-------------|------------|
| **Complexes** | 22 | 22 | ✓ Match |
| **Positive Examples** | 714 | 857 | +143 extra |
| **Negative Examples** | ~10,413-14,838 | 14,838 | ✓ Match |
| **Total** | ~11,127-15,552 | 15,695 | +143 |

**Key Finding**: The dataset file contains **143 more positive examples** than reported in the paper.

### Possible Explanations for Discrepancy

1. **Additional inhibitors** added after paper submission
2. **Filtered inhibitors** during quality control for the paper
3. **Duplicate removal** in paper analysis
4. **Complex filtering** - some complexes with only 1 inhibitor were removed

---

## Using the Preprocessed Dataset

### Step 1: Run the Preprocessing Script

```bash
python dataset_preprocessing_for_paper.py
```

This will:
- Load and analyze the main dataset
- Filter to the 22 complexes from the paper
- Separate external test sets
- Generate statistics and comparison with paper
- Save preprocessed files to `Data/preprocessed/`

### Step 2: Use the Preprocessed Files

Three files are generated:

#### 1. Training Dataset
**File**: `Data/preprocessed/training_22_complexes.txt`
- **Size**: 15,695 examples (857 positive, 14,838 negative)
- **Use**: Leave-One-Complex-Out (LOCO) cross-validation
- **Format**: Same as original file

#### 2. External Test Set 1 (2dyh - MDM2/p53)
**File**: `Data/preprocessed/external_test_2dyh.txt`
- **Size**: 72 examples (24 positive, 48 negative)
- **Use**: Independent validation on novel MDM2-p53 inhibitors

#### 3. External Test Set 2 (6m0j - SARS-CoV-2/ACE2)
**File**: `Data/preprocessed/external_test_6m0j.txt`
- **Size**: 72 examples (24 positive, 48 negative)
- **Use**: Independent validation on SARS-CoV-2 Spike/ACE2 inhibitors

---

## Feature Extraction Pipeline

### Overview

The paper uses a **three-stream feature extraction** approach:

```
┌─────────────────┐
│ Protein Complex │ ──→ GNN Features (512-dim)
│    (3D PDB)     │ ──→ Interface Features (211-dim)
└─────────────────┘ ──→ Sequence Features (69-dim)
                              ↓
┌─────────────────┐      Concatenate
│    Compound     │           ↓
│    (SMILES)     │ ──→ Morgan FP (2048-dim)
└─────────────────┘           ↓
                        Total: 2840-dim
                              ↓
                    ┌──────────────────┐
                    │ MLP Classifier   │
                    │ (512→100→1)      │
                    └──────────────────┘
```

### 1. Compound Features (2048-dim)

- **Method**: Extended-Connectivity Fingerprint (ECFP) / Morgan Fingerprint
- **Radius**: 2 bonds
- **Size**: 2048 bits
- **Tool**: RDKit

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def get_morgan_fingerprint(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    return np.array(fp)
```

### 2. Protein Complex Features (512-dim from GNN)

- **Method**: 3-layer Graph Neural Network (GNN)
- **Input**: 3D protein structure (PDB file)
- **Node Features**: Atom type (12-dim) + Residue type (21-dim)
- **Edge Definition**: Atoms within 6 Ångströms
- **Neighbor Aggregation**:
  - Same residue: 10 nearest neighbors
  - Different residue: 10 nearest neighbors
- **Output**: 512-dimensional embedding

### 3. Interface Features (211-dim)

- **Method**: Amino acid pair counts at interface
- **Interface Definition**: Residues within 8.0 Ångströms between chains
- **Encoding**: 20 × 20 amino acid pairs + 1 unknown = 211 dimensions
- **Normalization**: L2 normalization

```python
def extract_interface_features(pdb_file, distance_threshold=8.0):
    # Identify interface residues (distance < 8Å between chains)
    # Count amino acid pairs at interface
    # Return 211-dimensional vector
    pass
```

### 4. Protein Sequence Features (69-dim)

- **Amino Acid Composition (AAC)**: 20-dim (frequency of each amino acid)
- **Grouped 2-mer**: 49-dim (7 groups × 7 groups)
  - Groups based on physicochemical properties:
    - Aliphatic: {A, V, L, I, M, C}
    - Aromatic: {F, W, Y}
    - Polar-uncharged: {S, T, N, Q}
    - Basic: {K, R, H}
    - Acidic: {D, E}
    - Special: {G, P}
    - Other

---

## Leave-One-Complex-Out (LOCO) Validation

### Protocol

The paper uses **Leave-One-Complex-Out (LOCO)** cross-validation:

1. **For each of the 22 complexes**:
   - Remove all examples of that complex from training
   - Train model on remaining 21 complexes
   - Test on the held-out complex
   - Record AUC-ROC and AUC-PR

2. **Average results** across all 22 folds

3. **Per-complex results** reported in Table 2

### Expected Results (from Paper)

- **Average AUC-ROC**: 0.863 ± 0.096
- **Average AUC-PR**: 0.39 ± 0.236
- **Best complex**: 3DAB (AUC-ROC: 0.999)
- **Worst complex**: 2B4J (AUC-ROC: 0.663)

---

## Preprocessing Best Practices

### 1. Data Loading

```python
def load_dataset(filepath):
    """Load the PPI dataset."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            complex_name = parts[0]
            target_complex = parts[1]
            smiles = ' '.join(parts[2:-1])  # SMILES can have spaces
            label = float(parts[-1])

            data.append({
                'complex': complex_name,
                'target': target_complex,
                'smiles': smiles,
                'label': label
            })
    return data
```

### 2. Train/Test Split for LOCO

```python
def loco_split(data, test_complex):
    """Leave-One-Complex-Out split."""
    train = [x for x in data if x['complex'].split('_')[0] != test_complex]
    test = [x for x in data if x['complex'].split('_')[0] == test_complex]
    return train, test
```

### 3. Feature Extraction

```python
# Extract features for each example
for example in data:
    # 1. Compound features
    compound_fp = get_morgan_fingerprint(example['smiles'])

    # 2. Protein features (load precomputed or compute)
    protein_features = load_protein_features(example['complex'])

    # 3. Interface features (load precomputed or compute)
    interface_features = load_interface_features(example['complex'])

    # 4. Concatenate all features
    features = np.concatenate([
        compound_fp,           # 2048-dim
        protein_features,      # 512-dim
        interface_features     # 211-dim + 69-dim
    ])
```

### 4. Handle Class Imbalance

The paper uses **weighted loss** due to class imbalance:

```python
# Calculate class weights
pos_count = sum(1 for x in data if x['label'] == 1.0)
neg_count = len(data) - pos_count
pos_weight = neg_count / pos_count

# Use in loss function
criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([pos_weight]))
```

---

## Validation Checklist

Before training your model, verify:

- [ ] **22 complexes** are present in your training data
- [ ] **Positive examples** are properly labeled (label = 1.0)
- [ ] **Negative examples** are properly labeled (label = 0.0 or -1.0)
- [ ] **SMILES strings** are valid and can be parsed by RDKit
- [ ] **PDB files** exist for all protein complexes
- [ ] **Interface features** are precomputed or can be computed
- [ ] **Train/test split** ensures no complex overlap
- [ ] **Class weights** are computed for imbalanced data
- [ ] **External test sets** are kept separate and not used in training

---

## Troubleshooting

### Issue 1: Positive Example Count Mismatch

**Problem**: File has 857 positives, paper reports 714

**Solutions**:
1. **Use the full file** (857 positives) - More data may improve performance
2. **Filter to match paper** - Remove 143 examples based on:
   - Complexes with single inhibitors
   - Duplicate SMILES
   - Low-quality structures

### Issue 2: Missing PDB Files

**Problem**: Some protein complexes don't have PDB files in `Data/Pdb/`

**Solutions**:
1. Download from PDB: `https://www.rcsb.org/structure/{PDB_ID}`
2. Use precomputed features from Google Drive (see README)

### Issue 3: Invalid SMILES

**Problem**: Some SMILES strings cannot be parsed by RDKit

**Solutions**:
```python
from rdkit import Chem

def is_valid_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

# Filter invalid examples
data = [x for x in data if is_valid_smiles(x['smiles'])]
```

---

## References

1. **Main Paper**: Yaseen et al., "Predicting small-molecule inhibition of protein complexes", bioRxiv 2024
2. **2P2I Database**: Basse et al., "2P2Idb v2", Database 2016
3. **BindingDB**: Gilson et al., "BindingDB in 2015", Nucleic Acids Research 2016
4. **DBD5**: Vreven et al., "Updates to the Integrated Protein-Protein Interaction Benchmarks", JMB 2015
5. **SuperDRUG2**: Siramshetty et al., "SuperDRUG2", Nucleic Acids Research 2018

---

## Quick Start Commands

### 1. Preprocess the Dataset

```bash
cd /home/user/PPI-Inhibitors-main
python dataset_preprocessing_for_paper.py
```

### 2. Verify Preprocessing

```bash
# Check preprocessed files
ls -lh Data/preprocessed/

# View statistics
cat Data/preprocessed/dataset_statistics.txt
```

### 3. Load Preprocessed Data

```python
# In your notebook/script
from dataset_preprocessing_for_paper import PPIDatasetPreprocessor

preprocessor = PPIDatasetPreprocessor()
filtered_data, external_data, output_dir = preprocessor.run_full_preprocessing()

# Or load directly
train_data = pd.read_csv(
    'Data/preprocessed/training_22_complexes.txt',
    sep=' ',
    names=['complex', 'target', 'smiles', 'label']
)
```

---

## Summary

✅ **Use the preprocessed dataset** in `Data/preprocessed/training_22_complexes.txt`
✅ **22 complexes** match the paper
✅ **Leave-One-Complex-Out** validation protocol
✅ **External test sets** for independent validation
✅ **Three-stream features**: Compound (2048) + Protein GNN (512) + Interface+Sequence (280)
✅ **Weighted loss** for class imbalance

**Next Step**: Use the improved pipeline notebook to train the model with the preprocessed dataset!
