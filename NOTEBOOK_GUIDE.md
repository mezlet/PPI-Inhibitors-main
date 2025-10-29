# Comprehensive GNN-Based PPI Inhibitor Prediction Pipeline - Guide

## Overview

This document describes the comprehensive Jupyter notebook that replicates the complete GNN-based model from the research paper "Predicting small-molecule inhibition of protein complexes" by Yaseen et al. (2024).

## Notebook File

**Filename:** `Comprehensive_GNN_PPI_Inhibitor_Prediction_Pipeline.ipynb`

## Current Status

The notebook has been created with foundational sections including:

### ✅ Completed Sections

1. **Introduction and Overview** (Markdown)
   - Research paper summary
   - Problem statement and solution
   - Results overview
   - Notebook structure

2. **Section 1: Environment Setup** (Markdown + Code)
   - Mount Google Drive for pre-computed features
   - Clone GitHub repository
   - Install BioPython
   - Install RDKit

3. **Section 2: Import Libraries** (Code)
   - Deep learning libraries (PyTorch)
   - Structural biology (BioPython)
   - Cheminformatics (RDKit)
   - Machine learning (Scikit-learn)
   - Data manipulation (NumPy, Pandas)
   - CUDA configuration and utilities

4. **Section 3: Custom Data Loaders** (Markdown + Code)
   - Detailed explanation of class imbalance problem
   - `BalancedDataset` class with weighted sampling
   - `create_balanced_loader` function
   - `BinaryBalancedSampler` class for strict 50/50 splits
   - `CustomDataset` class
   - Testing and usage examples

5. **Section 4: Protein Feature Extraction** (Markdown - Header)
   - Overview of protein structure hierarchy
   - Explanation of features to be extracted

### 📝 Remaining Sections to Complete

The following sections need to be added to complete the notebook:

#### 6. **Protein Feature Extraction Functions (Code)**
```python
- atom1(structure): One-hot encode atom types (13 categories)
- res1(structure): One-hot encode residue types (21 amino acids)
- neigh1(structure): Build neighborhood graph (6Å threshold, 10 neighbors)
- GNN.processProtein(): Convert PDB files to GNN input format
```

#### 7. **Compound Feature Extraction**
```python
- Extract ECFP/Morgan fingerprints from SMILES
- Radius=2, 2048-bit fingerprint
- Using RDKit
```

#### 8. **Sequence and Interface Features**
```python
- Amino acid composition (AAC): 20-dimensional
- Grouped 2-mer features: 49-dimensional (7 physicochemical groups)
- Interface residue pairs: 211-dimensional (residues within 8Å)
```

#### 9. **GNN Model Architecture**
```python
- GNN_First_Layer: Initial layer with atom and residue encoding
- GNN_Layer: Intermediate graph convolution layers
- Dense: Output layer
- GNN: Complete model (3 conv layers: 512 → 1024 → 512)
```

#### 10. **MLP Classifier**
```python
- IPPI_MLP_Net: Multi-layer perceptron
- Input: 2840 dimensions (512 GNN + 280 interface/sequence + 2048 compound)
- Architecture: 2840 → 1024 → 512 → 100 → 1
- Activations: tanh, tanh, ReLU, none
```

#### 11. **Data Loading and Preparation**
```python
- Load 2P2I positive examples (714 inhibitors)
- Load negative examples (10,413 total):
  * Random 2P2I/SuperDRUG2 pairs
  * 2P2I compounds with DBD5 complexes
  * Binders with Tanimoto < 0.85
- Load pre-computed features
- Create data dictionaries
```

#### 12. **Training Pipeline - Leave-One-Complex-Out CV**
```python
- GroupKFold with complex names as groups
- For each fold:
  * Hold out one complex as test set
  * Train on remaining 21 complexes
  * Standardize features
  * Create balanced data loaders
  * Initialize GNN and MLP models
  * Training loop with early stopping
  * Save best model based on validation AUC-ROC
```

#### 13. **Training Configuration**
```python
- Loss: Binary Cross Entropy with Logits
- Optimizer: Adam (lr=0.0001, weight_decay=0.0)
- Batch size: 1024
- Weighted loss (inverse class frequencies per complex)
- Early stopping based on validation AUC-ROC
- GPU acceleration
```

#### 14. **Model Evaluation**
```python
- Calculate metrics for each held-out complex:
  * AUC-ROC
  * AUC-PR
  * Precision-Recall curves
  * ROC curves
- Aggregate results across all folds
- Report mean ± std for each complex
```

#### 15. **External Test Set 1: Recent Publications**
```python
- Load 28 inhibitors from recent literature
- Generate negative examples
- Test best models from CV
- Calculate AUC-ROC and AUC-PR
- Expected result: AUC-ROC ~0.82
```

#### 16. **External Test Set 2: SARS-CoV-2 Inhibitors**
```python
- Load 25 RBD-hACE2 PPI inhibitors
- Generate negative examples
- Test best models from CV
- Calculate AUC-ROC and AUC-PR
- Expected result: AUC-ROC ~0.78
```

#### 17. **Baseline Comparison: SVM**
```python
- Load SVM predictions
- Calculate SVM metrics
- Expected: AUC-ROC ~0.74, AUC-PR ~0.33
```

#### 18. **Baseline Comparison: GearNet Embeddings**
```python
- Load GearNet embedding predictions
- Calculate GearNet metrics
- Expected: AUC-ROC ~0.79, AUC-PR ~0.30
```

#### 19. **Results Visualization**
```python
- Combined ROC curves (GNN, SVM, GearNet)
- Combined PR curves (GNN, SVM, GearNet)
- Per-complex performance barplots
- Save figures as PDF
```

#### 20. **Results Summary and Discussion**
```markdown
- Table of per-complex results
- Overall performance comparison
- Discussion of findings
- Limitations and future work
```

## Key Features of This Notebook

### 1. **Comprehensive Documentation**
- Every function has detailed docstrings
- Mathematical formulations explained
- Code comments on nearly every line
- Markdown cells explain concepts before code

### 2. **Modular Structure**
- Clear separation of concerns
- Reusable functions
- No code duplication
- Follows original research methodology

### 3. **Educational Value**
- Explains the "why" not just the "what"
- Examples and usage demonstrations
- Links to research paper concepts
- Suitable for learning and reproduction

### 4. **Complete Reproduction**
- Uses exact functions from original codebase
- Same hyperparameters
- Same data preprocessing
- Same evaluation methodology

## Data Requirements

### Pre-computed Features (Google Drive)
These files must be downloaded from the Google Drive links in the README:

1. **ProteinData_dict.pickle** (~1-2 GB)
   - Pre-computed GNN features for 2P2I complexes
   - Link: https://drive.google.com/file/d/1goeDiPZSKT1Xx3j00eNG9xlqYkLLv1gW/view

2. **DBD5_ProteinData_dict.pickle** (~1-2 GB)
   - Pre-computed GNN features for DBD5 complexes
   - Link: https://drive.google.com/file/d/1GOYEKLQCoGea9QQ72kujy0rdJKbUSYAE/view

### Repository Data (Included)
These are available in the cloned repository:

1. **2p2iComplexPairs.txt**
   - Protein complex sequences

2. **2p2iInhibitorsSMILES.txt**
   - Inhibitor SMILES strings

3. **WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt**
   - Complete dataset with labels

4. **Binders With Tanimoto Similarity 0.85.csv**
   - Hard negative examples (binders that are not inhibitors)

5. **approved_drugs_chemical_structure_identifiers.xlsx**
   - SuperDRUG2 compounds for random negatives

### Pre-computed Feature Dictionaries
Located in `Features/` directory:

1. **Compound_Fingerprint_Features_Dict.npy**
   - ECFP fingerprints for all compounds

2. **Pos_seqandInterfaceF_dict.npy**
   - Sequence and interface features for 2P2I

3. **NewUbench5InterfaceandSeq_dict.npy**
   - Sequence and interface features for DBD5

4. **Classratio_GNNdict.npy**
   - Class ratios for weighted loss

## Expected Runtime

On Google Colab with GPU:
- Setup and data loading: ~5-10 minutes
- Training (22-fold CV): ~4-6 hours
- External testing: ~10-20 minutes
- **Total**: ~5-7 hours

On CPU only:
- **Total**: ~24-48 hours (not recommended)

## Expected Results

### Cross-Validation (LOCO)
- Mean AUC-ROC: 0.863 ± 0.096
- Mean AUC-PR: 0.39 ± 0.236

### External Test 1 (Recent Publications)
- AUC-ROC: 0.82

### External Test 2 (COVID-19)
- AUC-ROC: 0.76-0.78

### Baseline Comparisons
| Method | AUC-ROC | AUC-PR |
|--------|---------|--------|
| **GNN (Proposed)** | **0.863** | **0.39** |
| SVM Kernel | 0.744 | 0.33 |
| GearNet Embeddings | 0.794 | 0.30 |

## How to Use This Notebook

### On Google Colab (Recommended)

1. **Open the notebook in Colab**
   - Click the "Open in Colab" badge at the top

2. **Set up GPU runtime**
   - Runtime → Change runtime type → GPU → Save

3. **Mount Google Drive**
   - Run first cell
   - Authenticate your Google account
   - Download pre-computed features to: `/content/drive/MyDrive/GNN-PPI-Inhibitor/`

4. **Run all cells sequentially**
   - Runtime → Run all
   - Or execute cells one by one to follow along

5. **Monitor progress**
   - Training progress shown with tqdm bars
   - Validation results printed after each complex

### On Local Machine

1. **Prerequisites**
   - Python 3.6+
   - NVIDIA GPU with CUDA (highly recommended)
   - 16+ GB RAM
   - 10+ GB disk space

2. **Install dependencies**
   ```bash
   pip install torch torchvision biopython rdkit scikit-learn pandas numpy matplotlib tqdm
   ```

3. **Download data**
   - Clone repository
   - Download Google Drive files
   - Place in appropriate directories

4. **Update paths**
   - Change `path` and `githubpath` variables
   - Point to your local directories

5. **Run notebook**
   ```bash
   jupyter notebook Comprehensive_GNN_PPI_Inhibitor_Prediction_Pipeline.ipynb
   ```

## Troubleshooting

### Common Issues

**1. CUDA Out of Memory**
```
Solution: Reduce batch_size from 1024 to 512 or 256
```

**2. Missing Features**
```
Error: FileNotFoundError: ProteinData_dict.pickle
Solution: Download from Google Drive links provided
```

**3. PDB Warnings**
```
Warning: PDBConstructionWarning
Solution: These are normal, the code catches them
```

**4. Slow Training**
```
Issue: Training is very slow
Solution: Ensure GPU is enabled in Colab runtime settings
```

## Code Style and Conventions

### Variable Naming
- `P`: Protein features
- `C`: Compound features
- `Y`: Labels
- `Z`: Predictions
- `tr`: Training set
- `tt`: Test set
- `dict`: Dictionary data structure

### Function Naming
- `atom1()`, `res1()`, `neigh1()`: Feature extraction
- `processProtein()`: Batch processing
- `cuda()`, `toTensor()`, `toNumpy()`: Type conversions

### File Naming
- Models saved as: `GNN-based-pipeline_IPPI_Net_{complex_name}`
- Results saved as: `{complex_name}Scores.npy`, `{complex_name}Targets.npy`

## Citation

If you use this notebook or the methodology, please cite the original paper:

```bibtex
@article{yaseen2024predicting,
  title={Predicting small-molecule inhibition of protein complexes},
  author={Yaseen, Adiba and Roy, Soumyadip and Akhter, Naeem and Ben-Hur, Asa and Minhas, Fayyaz},
  journal={bioRxiv},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Repository Structure

```
PPI-Inhibitors-main/
├── Comprehensive_GNN_PPI_Inhibitor_Prediction_Pipeline.ipynb  # THIS NOTEBOOK
├── NOTEBOOK_GUIDE.md                                          # THIS GUIDE
├── README.md                                                  # Original README
├── research_paper.pdf                                         # Research paper
├── Data/
│   ├── 2p2iComplexPairs.txt
│   ├── 2p2iInhibitorsSMILES.txt
│   ├── WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt
│   ├── Binders With Tanimoto Similarity 0.85.csv
│   └── approved_drugs_chemical_structure_identifiers.xlsx
├── Features/
│   ├── Compound_Fingerprint_Features_Dict.npy
│   ├── Pos_seqandInterfaceF_dict.npy
│   ├── NewUbench5InterfaceandSeq_dict.npy
│   └── Classratio_GNNdict.npy
├── code/
│   ├── GNN_based_pipeline_Training_for_Predicting...ipynb     # Original training
│   └── seqfeaturesand_gnn_generate_prediction...ipynb         # Original testing
└── Final Results/
    └── Figures/
```

## Next Steps for Completion

To complete the notebook, the following cells need to be added:

1. Add protein feature extraction code (Section 6)
2. Add compound feature extraction code (Section 7)
3. Add sequence/interface features code (Section 8)
4. Add GNN model definition (Section 9)
5. Add MLP model definition (Section 10)
6. Add data loading code (Section 11)
7. Add training loop (Section 12-13)
8. Add evaluation code (Section 14)
9. Add external testing (Sections 15-16)
10. Add baseline comparisons (Sections 17-18)
11. Add visualization code (Section 19)
12. Add final summary (Section 20)

Each section should follow the same documentation style as Sections 1-4:
- Detailed markdown explanations before code
- Comprehensive code comments
- Example usage where appropriate
- Mathematical formulations explained
- Links to research paper concepts

## Contact and Support

For questions about this notebook:
- Open an issue on the GitHub repository
- Refer to the original research paper for methodology details
- Check the README.md for dataset information

---

**Document Version:** 1.0
**Last Updated:** October 29, 2025
**Author:** Based on work by Yaseen et al. (2024)
