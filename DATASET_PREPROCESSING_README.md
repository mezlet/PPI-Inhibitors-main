# PPI Inhibitor Dataset Preprocessing Guide

## Overview

This document describes how to preprocess the PPI inhibitor dataset to exactly match the specifications from the research paper:

**"Predicting small-molecule inhibition of protein complexes"** by Yaseen et al.

## Research Paper Dataset Specifications

### Positive Examples (Section 2.1.1)

**Source:** 2P2I v2 database

**Filtering Pipeline:**
1. **Initial dataset:** 822 protein complex-inhibitor pairs from 32 complexes
2. **Remove predicted structures:** 7 complexes with only predicted (not resolved) structures removed → **722 examples from 25 complexes**
3. **Remove single-inhibitor complexes:** Complexes with only 1 inhibitor removed (for robust performance assessment) → **714 examples**
4. **Final result:** **714 inhibitors across 22 complexes with 608 unique inhibitor compounds**

### Negative Examples (Section 2.1.2)

Three strategies producing **10,413 negative examples** total:

#### 1. Random 2P2I + SuperDRUG2 Pairing (~857 examples)
- 2P2I complexes paired with random compounds from:
  - 2P2I database
  - SuperDRUG2 database (FDA-approved drugs)
- Total pool: 3,633 unique small molecules
- **Constraint:** Selected compound must NOT be a known inhibitor of that complex

#### 2. 2P2I Compounds + DBD5 Complexes (~1,714 examples)
- 2P2I inhibitor compounds paired with DBD5 complexes
- **DBD5:** Protein-protein docking benchmark database v5.5
- Total: 282 complexes from DBD5
- **Constraint:** Only complexes with bound 3D structures used

#### 3. Binders that are NOT Inhibitors (~7,842 examples)
- **Source:** BindingDB database
- **Strategy:** Find compounds that bind protein chains but do NOT inhibit the complex
- **Selection Criteria:**
  1. **BLASTp search:** >90% protein sequence identity with 2P2I chains
  2. **Strong binding affinity:** Ki, Kd, or IC50 < 7.6 nM (active binders)
  3. **Tanimoto dissimilarity:** <0.85 with any known inhibitor of that complex (excludes possible inhibitors)
- **Purpose:** Creates "hard" negative examples (binders that don't inhibit)

### Final Dataset Summary

| Category | Count |
|----------|-------|
| **Positive Examples** | 714 |
| **Negative Examples** | 10,413 |
| **Total Examples** | 11,127 |
| **Unique Complexes** | 22 |
| **Unique Compounds (positive)** | 608 |
| **Positive:Negative Ratio** | 1:14.6 |

---

## Current Dataset Files

### Main Dataset File
**`WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt`**

**Format:**
```
[TestPositiveComplex] [TargetComplex] [CompoundName] [Label]
```

**Example:**
```
3UVW_A_2_B 3UVW_A_2_B WSH 1.0
3UVW_A_2_B 1BXL_A_2_B ASPIRIN 0.0
```

**Current Statistics:**
- Total examples: 15,695
- Positive (label=1.0): 857
- Negative (label=0.0): 14,838
- Unique complexes (positive): 22

**Note:** This file contains MORE positive examples (857) than the paper's final dataset (714), which suggests it includes all inhibitors before the final filtering step described in the paper.

---

## Preprocessing Pipeline

### Notebook: `Dataset_Preprocessing_Pipeline.ipynb`

This notebook filters the raw dataset to match the paper's exact specifications.

### Preprocessing Steps

```
Raw Dataset (15,695 examples, 857 positive, 22 complexes)
    ↓
Step 1: Remove complexes with only 1 inhibitor
    ↓
Step 2: Keep top 22 complexes (matches paper's complex count)
    ↓
Step 3: Sample exactly 714 positive examples (matches paper)
    ↓
Step 4: Keep all negative examples
    ↓
Filtered Dataset (matching paper specifications)
```

### Output Files

1. **`WriteAllexamples_Filtered_Paper_Specs.txt`**
   - Filtered dataset matching paper specifications
   - Format: same as input file
   - Ready for use with GNN pipeline

2. **`Filtered_Complex_List.txt`**
   - List of 22 complexes in filtered dataset
   - Format: `ComplexID  NumInhibitors`

3. **`Dataset_Statistics.txt`**
   - Detailed statistics comparing raw, filtered, and paper datasets

4. **`Preprocessing_Report.txt`**
   - Comprehensive preprocessing report

---

## Usage Instructions

### 1. Run Preprocessing Notebook

```bash
# If using Google Colab:
1. Upload Dataset_Preprocessing_Pipeline.ipynb to Colab
2. Clone the repository or upload data files
3. Run all cells
4. Download generated files

# If using local Jupyter:
jupyter notebook Dataset_Preprocessing_Pipeline.ipynb
```

### 2. Update Main Pipeline

In `Complete_PPI_Inhibitors_Pipeline_End_To_End.ipynb`, update the data loading cell:

```python
# OLD:
input_file = githubpath + 'Data/WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt'

# NEW:
input_file = githubpath + 'Data/WriteAllexamples_Filtered_Paper_Specs.txt'
```

### 3. Run Training

Execute the main pipeline with the filtered dataset. Expected results should match the paper:

**Cross-Validation (Leave-One-Complex-Out):**
- AUC-ROC: 0.86 ± 0.10
- AUC-PR: 0.39 ± 0.24

**External Validation:**
- Recent publications dataset: AUC-ROC 0.82
- SARS-CoV-2 (RBD-hACE2): AUC-ROC 0.78

---

## Dataset Construction Details (from Paper)

### Feature Extraction

#### 1. Compound Features (2048-dimensional)
- **Method:** Extended-Connectivity Fingerprint (ECFP) / Morgan Fingerprint
- **Tool:** RDKit
- **Parameters:**
  - Radius: 2 bonds
  - Number of bits: 2048
- **Input:** SMILES representation

#### 2. Protein Sequence Features (69-dimensional)
- **Amino Acid Composition (AAC):** 20-dimensional
  - Frequency of each amino acid
- **Grouped k-mer (k=2):** 49-dimensional
  - Amino acids grouped by physicochemical properties (7 groups)
  - Counts of grouped k-mers
  - 7² = 49 features

#### 3. Interface Features (211-dimensional)
- **Definition:** Residues within 8Å of each other across chains
- **Features:** Frequency of amino acid pairs at interface
- **Dimension:** 21 × 21 = 441, but condensed to 211 unique pairs

#### 4. GNN Features (512-dimensional)
- **Input:** 3D protein complex structure
- **Architecture:** 3-layer heterogeneous GNN
  - Layer 1: 512 dimensions (atom + residue features)
  - Layer 2: 1024 dimensions
  - Layer 3: 512 dimensions
- **Output:** Global pooled graph embedding

#### 5. Combined Features (2840-dimensional)
- GNN features: 512
- Interface features: 211
- Sequence features: 69
- Compound features: 2048
- **Total:** 2840 dimensions → passed to MLP

### Model Architecture

```
Input: (Protein Complex, Compound)
    ↓
┌─────────────────┬──────────────────┐
│   GNN Branch    │  Feature Branch  │
│   (3D structure)│  (pre-computed)  │
│                 │                  │
│  Atom features  │  Interface (211) │
│  Residue feats  │  Sequence (69)   │
│  Neighbors      │  Compound (2048) │
│       ↓         │        ↓         │
│  GNN Layer 1    │   Concatenate    │
│  (512 dims)     │                  │
│       ↓         │                  │
│  GNN Layer 2    │                  │
│  (1024 dims)    │                  │
│       ↓         │                  │
│  GNN Layer 3    │                  │
│  (512 dims)     │                  │
│       ↓         │                  │
│  Global Pool    │                  │
│  (512 dims)     │                  │
└─────────────────┴──────────────────┘
            ↓
    Concatenate (2840 dims)
            ↓
    MLP Layer 1 (1024, tanh)
            ↓
    MLP Layer 2 (512, tanh)
            ↓
    MLP Layer 3 (100, ReLU)
            ↓
    Output Layer (1, sigmoid)
            ↓
    Inhibition Score [0, 1]
```

### Training Strategy

```python
# Hyperparameters (from paper)
optimizer = Adam(lr=0.0001, weight_decay=0.0)
loss_function = BCEWithLogitsLoss()
batch_size = 1024
epochs = 2  # (adjustable)

# Balanced Sampling
# - Each batch: 50% positive, 50% negative
# - Undersample majority class to match minority

# Weighted Loss
# - Weight by positive:negative ratio per complex
# - Accounts for class imbalance

# Validation Strategy
# - Leave-One-Complex-Out (LOCO) Cross-Validation
# - 22 folds (one per complex)
# - Train on 21 complexes, test on 1 held-out complex
```

---

## Validation Protocols (from Paper)

### 1. Cross-Validation (LOCO)

**Purpose:** Assess generalization to unseen protein complexes

**Method:**
- **Leave-One-Complex-Out (LOCO) Cross-Validation**
- **Folds:** 22 (one per complex)
- **Training:** All examples from 21 complexes
- **Testing:** All examples from 1 held-out complex

**Rationale:**
- Tests ability to predict inhibitors for completely novel protein complexes
- More rigorous than random splits (avoids data leakage)
- Mimics real-world scenario: predicting inhibitors for new targets

**Expected Results:**
- Mean AUC-ROC: 0.86 (SD: 0.10)
- Mean AUC-PR: 0.39 (SD: 0.24)

### 2. External Validation

#### Dataset 1: Recent Publications
- **Source:** 28 inhibitors from recent literature (post-2P2I)
- **Novelty:** Different structures from training data
- **Result:** AUC-ROC 0.82

#### Dataset 2: SARS-CoV-2 Inhibitors
- **Source:** 25 inhibitors of RBD-hACE2 PPI (Hanson et al. 2020)
- **Complex:** SARS-CoV-2 Spike protein RBD + Human ACE2
- **Result:** AUC-ROC 0.78

---

## Key Differences: Raw vs Filtered Dataset

| Feature | Raw Dataset | Filtered Dataset | Paper Specs |
|---------|-------------|------------------|-------------|
| **Total Examples** | 15,695 | ~11,127 | 11,127 |
| **Positive** | 857 | 714 | 714 |
| **Negative** | 14,838 | ~10,413 | 10,413 |
| **Complexes** | 22 | 22 | 22 |
| **Unique Compounds** | ~800 | 608 | 608 |
| **Pos:Neg Ratio** | 1:17.3 | 1:14.6 | 1:14.6 |

**Why the difference?**
- The raw dataset includes ALL inhibitors found in 2P2I for the 22 complexes (857 total)
- The paper filtered further to remove:
  - Examples with predicted structures
  - Duplicate or low-quality examples
- Final result: 714 examples matching paper specifications

---

## Troubleshooting

### Issue: Results don't match paper

**Possible causes:**
1. **Wrong dataset file:** Ensure you're using the filtered dataset
2. **Different random seed:** Results may vary slightly due to random sampling
3. **Hyperparameters:** Verify all hyperparameters match the paper
4. **Pre-computed features:** Ensure you're using the correct pre-computed feature files

**Solutions:**
- Run the preprocessing notebook to regenerate the filtered dataset
- Check that the filtered dataset has exactly 714 positive examples and 22 complexes
- Verify that the feature files match the paper's specifications

### Issue: Missing pre-computed features

**Required feature files** (should be in `Features/` directory):
1. `Pos_seqandInterfaceF_dict.npy` - Interface features for 2P2I complexes
2. `NewUbench5InterfaceandSeq_dict.npy` - Interface features for DBD5 complexes
3. `Compound_Fingerprint_Features_Dict.npy` - Morgan fingerprints for all compounds
4. `Classratio_GNNdict.npy` - Class ratios for balanced sampling

**If missing:**
- Download from Google Drive links in README.md
- Or re-compute using the feature extraction code in the main pipeline

### Issue: Memory errors during training

**Solutions:**
1. **Reduce batch size:** Try batch_size=512 or 256
2. **Use GPU:** Training on GPU significantly reduces memory usage
3. **Clear cache:** Add `torch.cuda.empty_cache()` between training iterations

---

## References

**Paper:**
Yaseen, A., Roy, S., Akhter, N., Ben-Hur, A., & Minhas, F. (2024).
"Predicting small-molecule inhibition of protein complexes."
*Bioinformatics*, bioRxiv preprint doi: 10.1101/2024.08.23.609286

**Databases:**
- **2P2I v2:** Basse et al. (2016) - Protein-Protein Interaction inhibitor database
- **DBD5:** Vreven et al. (2015) - Protein docking benchmark database
- **BindingDB:** Gilson et al. (2016) - Binding affinity database
- **SuperDRUG2:** Siramshetty et al. (2018) - FDA-approved drugs database

**Code Repository:**
https://github.com/adibayaseen/PPI-Inhibitors

---

## Questions?

For questions about the preprocessing pipeline or dataset specifications:
1. Check the paper's Methods section (Section 2)
2. Review the Supplementary Materials
3. Examine the preprocessing notebook comments
4. Open an issue on the GitHub repository

---

## License

This preprocessing pipeline is part of the PPI-Inhibitors project and follows the same license as the original repository.
