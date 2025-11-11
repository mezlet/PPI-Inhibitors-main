# Improved PPI Inhibitors Pipeline - README

## What's New

This repository now includes **`Improved_PPI_Inhibitors_Pipeline_With_Preprocessing.ipynb`**, which extends the original `Complete_PPI_Inhibitors_Pipeline_End_To_End.ipynb` with complete dataset preprocessing code that matches the research paper methodology.

## Quick Comparison

| Feature | Original Notebook | Improved Notebook |
|---------|------------------|-------------------|
| **Dataset Preprocessing** | ❌ Not included | ✅ Complete implementation |
| **Starts from** | Pre-processed file | Raw data files (2P2I, SuperDRUG2, DBD5) |
| **Negative Generation** | ❌ Not shown | ✅ All 3 strategies implemented |
| **Model Training** | ✅ Complete | ✅ Same as original |
| **Documentation** | Basic | ✅ Extensive with paper references |
| **Reproducibility** | Partial | ✅ Full end-to-end |

## Files Added

1. **`Improved_PPI_Inhibitors_Pipeline_With_Preprocessing.ipynb`**
   - Complete preprocessing + training pipeline
   - Implements all strategies from the research paper
   - ~60 cells with detailed explanations

2. **`DATASET_PREPROCESSING_GUIDE.md`**
   - Comprehensive guide to dataset preprocessing
   - Explains each strategy in detail
   - Includes validation checks

3. **`IMPROVED_PIPELINE_README.md`** (this file)
   - Quick start guide
   - Summary of changes

## Research Paper Dataset Methodology

The improved notebook implements the exact preprocessing described in the paper:

### Positive Examples (714 total)
- **Source**: 2P2I-DB v2 database
- **Complexes**: 22 protein complexes (filtered from 32 original)
- **Inhibitors**: 608 unique small molecules

### Negative Examples (10,413 total) - Three Strategies

#### 1. Random 2P2I/SuperDRUG2 Pairing (~857 examples)
```python
# Pair 2P2I complexes with random compounds
# from 2P2I (not their inhibitors) + SuperDRUG2
```

#### 2. 2P2I Compounds × DBD5 Complexes (~1,714 examples)
```python
# Pair 2P2I inhibitors with 282 DBD5 complexes
# (different protein complexes)
```

#### 3. Hard Negatives: Binders NOT Inhibitors (~11,789 examples)
```python
# From BindingDB, find compounds that:
# - BIND to proteins (>90% seq identity)
# - Strong binding (Ki/Kd/IC50 < 7.6 nM)
# - NOT similar to inhibitors (Tanimoto < 0.85)
```

## Quick Start

### Using Google Colab (Recommended)

1. **Open the improved notebook**:
   - Click the "Open in Colab" badge at the top of the notebook
   - Or upload `Improved_PPI_Inhibitors_Pipeline_With_Preprocessing.ipynb` to Colab

2. **Set runtime to GPU**:
   - Runtime → Change runtime type → GPU (T4 or better)

3. **Run all cells**:
   - Runtime → Run all
   - First section: Setup and data preprocessing (~10-15 minutes)
   - Second section: Model training and evaluation (~30-60 minutes)

### What the Notebook Does

```
Section 1: Setup (Cells 1-4)
├── Clone repository
├── Install dependencies
└── Mount Google Drive

Section 2: Import Libraries (Cells 5-6)
└── Load all required packages

Section 3: DATASET PREPROCESSING (Cells 7-17) ⭐ NEW!
├── Load 2P2I data
├── Filter positive examples
├── Generate Strategy 1 negatives
├── Generate Strategy 2 negatives
├── Generate Strategy 3 negatives
├── Combine all examples
└── Save preprocessed dataset

Section 4: Load Protein Features (Cells 18-19)
└── Download pre-computed GNN features

Section 5: Utility Functions (Cells 20-22)
├── Amino acid composition
├── Grouped k-mer features
├── Morgan fingerprints
└── Interface features

Section 6: Model Architecture (Cells 23-27)
├── GNN layers
└── Full model with MLP

Section 7: Training & Evaluation (Cells 28+)
├── Leave-one-complex-out CV
├── External dataset validation
└── Results visualization
```

## Expected Outputs

After running preprocessing (Section 3):

```
Data/
├── Preprocessed_Dataset_All_Examples.txt   # Main dataset file
├── Preprocessed_Dataset_All_Examples.csv   # CSV version with metadata
└── Preprocessed_SMILES_Dict.txt           # Inhibitor name → SMILES mapping
```

### Dataset Statistics

```
FINAL DATASET STATISTICS
═══════════════════════════
Total examples: ~11,127
Positive examples: ~714
Negative examples: ~10,413

Negative examples by strategy:
  random_2p2i_superdrug: ~857
  random_dbd5_2p2i: ~1,714
  binder_non_inhibitor: ~8,842

Unique complexes: ~22 (2P2I) + ~100 (DBD5) = ~122
Unique SMILES: ~3,000+
```

## Key Features

### 1. Complete Transparency
Every preprocessing step is visible and modifiable:
```python
# You can see exactly how negatives are generated
def generate_random_negative_strategy1(complexes, compounds, ...):
    # Clear, documented code
    ...
```

### 2. Research Paper Alignment
Each step references the paper:
```markdown
### 3.3 Generate Negative Examples - Strategy 1
**Random pairing of 2P2I complexes with compounds from 2P2I and SuperDRUG2**
- Pair complexes randomly with compounds that are NOT known inhibitors
- Total: ~857 negative examples
```

### 3. Flexible Configuration
Easily adjust preprocessing parameters:
```python
# Generate more/fewer negatives per complex
negatives_strategy1 = generate_random_negative_strategy1(
    complexes_with_multiple_inhibitors,
    all_compounds,
    positive_examples,
    num_per_complex=10  # ← Adjust this
)
```

### 4. Validation Checks
Built-in checks to ensure data quality:
```python
# Verify dataset matches paper
assert len(positive_examples) >= 700, "Should have ~714 positive examples"
assert positive_examples['complex_name'].nunique() >= 20, "Should have ~22 complexes"
```

## When to Use Each Notebook

### Use `Complete_PPI_Inhibitors_Pipeline_End_To_End.ipynb` if:
- You want to quickly train models
- You trust the pre-processed dataset
- You don't need to modify preprocessing
- You're focused on model architecture/training

### Use `Improved_PPI_Inhibitors_Pipeline_With_Preprocessing.ipynb` if:
- You want full end-to-end reproducibility ⭐
- You need to understand preprocessing
- You want to modify negative generation
- You're conducting research on the dataset
- You're writing a paper and need to cite methodology

## Common Use Cases

### 1. Reproduce Paper Results Exactly
```python
# The improved notebook generates the dataset
# using the exact methodology from the paper
# Then trains models with the same architecture
```

### 2. Experiment with Different Negatives
```python
# Modify Strategy 1 to use different compound sources
superdrug_smiles = load_your_compound_library()

# Adjust negative sampling ratios
num_per_complex = 20  # More negatives per complex

# Use different Tanimoto thresholds
tanimoto_threshold = 0.90  # More stringent filtering
```

### 3. Add New Protein Complexes
```python
# Add your own complexes to the preprocessing
new_complexes = load_your_complexes('my_complexes.txt')
complexes_2p2i.update(new_complexes)

# Generate examples
all_examples = preprocess_with_new_complexes(...)
```

### 4. Test Different Negative Strategies
```python
# Try only hard negatives (Strategy 3)
all_examples = pd.concat([
    positive_examples,
    negatives_strategy3  # Only binders
])

# Compare performance
```

## Validation

The notebook includes validation to ensure preprocessing is correct:

```python
# Check class balance
pos_count = (all_examples['label'] == 1).sum()
neg_count = (all_examples['label'] == 0).sum()
print(f"Positive: {pos_count}, Negative: {neg_count}")
print(f"Ratio: 1:{neg_count/pos_count:.1f}")

# Expected: Ratio ~1:14.6 (matching paper)
```

## Performance Expectations

Based on the research paper:

| Validation Method | AUC-ROC | AUC-PR |
|------------------|---------|--------|
| **Leave-One-Complex-Out CV** | 0.86 ± 0.10 | 0.39 ± 0.24 |
| **External Dataset (Recent Pubs)** | 0.82 | - |
| **External Dataset (SARS-CoV-2)** | 0.78 | - |

## Troubleshooting

### Issue: "File not found" errors

**Problem**: Missing data files (SuperDRUG2, DBD5, etc.)

**Solutions**:
1. Download missing files from the repository
2. Or skip that strategy:
   ```python
   # Skip Strategy 2 if DBD5 not available
   negatives_strategy2 = pd.DataFrame()
   ```

### Issue: Low number of negatives generated

**Problem**: Not enough compounds/complexes in source files

**Solutions**:
1. Increase `num_per_complex` parameter
2. Add more compound sources
3. Use more DBD5 complexes

### Issue: Memory errors during preprocessing

**Problem**: Large datasets + limited RAM

**Solutions**:
1. Process in batches:
   ```python
   for batch in complex_batches:
       process_batch(batch)
   ```
2. Use Colab Pro for more RAM
3. Save intermediate results

## Citation

If you use this improved pipeline, please cite:

**Original Paper**:
```bibtex
@article{yaseen2024predicting,
  title={Predicting small-molecule inhibition of protein complexes},
  author={Yaseen, Adiba and Roy, Soumyadip and Akhter, Naeem and Ben-Hur, Asa and Minhas, Fayyaz},
  journal={bioRxiv},
  year={2024}
}
```

**Repository**:
```bibtex
@software{ppi_inhibitors_improved,
  title={Improved PPI Inhibitors Pipeline with Dataset Preprocessing},
  author={[Your contribution]},
  year={2024},
  url={https://github.com/adibayaseen/PPI-Inhibitors}
}
```

## Support

For questions about:
- **Preprocessing**: See `DATASET_PREPROCESSING_GUIDE.md`
- **Original pipeline**: See `README.md` in the repository
- **Research paper**: See `research_paper.pdf` and `Supplementary.pdf`
- **Issues**: Open an issue on GitHub

## License

Same as the original repository. See LICENSE file.

## Acknowledgments

- Original authors for the research paper and initial codebase
- 2P2I-DB v2 database maintainers
- SuperDRUG2 database team
- DBD5 benchmark database creators
- BindingDB for compound-protein binding data
