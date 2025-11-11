# PPI Inhibitors Dataset Preprocessing Guide

## Overview

This guide explains how the dataset preprocessing in `Improved_PPI_Inhibitors_Pipeline_With_Preprocessing.ipynb` matches the methodology described in the research paper: **"Predicting small-molecule inhibition of protein complexes"**.

## Research Paper Dataset Description

According to the research paper (Section 2.1), the dataset consists of:

### Positive Examples (714 examples from 22 complexes)
- **Source**: 2P2I-DB v2 database
- **Initial size**: 822 inhibitors from 32 complexes
- **Filtering applied**:
  - Removed 7 complexes with only predicted structures (not experimentally verified)
  - Removed complexes with only 1 inhibitor
- **Final dataset**: 714 positive examples from 22 protein complexes with 608 unique inhibitors

### Negative Examples (10,413 total examples)

The paper uses three strategies to generate negative examples:

#### Strategy 1: Random Pairing with 2P2I and SuperDRUG2 (~857 examples)
- Randomly pair complexes from 2P2I with compounds from:
  - Other 2P2I inhibitors (not inhibitors of that complex)
  - SuperDRUG2 approved drugs database (3,633 unique compounds)
- Rationale: Random compound-complex pairing is unlikely to be inhibitory

#### Strategy 2: 2P2I Compounds with DBD5 Complexes (~1,714 examples)
- Pair 2P2I inhibitors with 282 complexes from DBD5 benchmark database (v5.5)
- Use only complexes with bound 3D structures available
- Rationale: Compounds from one complex are unlikely to inhibit unrelated complexes

#### Strategy 3: Binders that are NOT Inhibitors (~11,789 examples)
This is the most sophisticated "hard negative" strategy:

1. **BLASTp search in BindingDB**:
   - For each chain in each complex
   - Find proteins with >90% sequence identity

2. **Filter for strong binders**:
   - Keep only ligands with measured binding affinity
   - Criteria: Ki, Kd, or IC50 < 7.6 nM (strong binders)

3. **Exclude similar inhibitors**:
   - Calculate Tanimoto coefficient with known inhibitors
   - Keep only binders with Tanimoto < 0.85 (dissimilar to known inhibitors)

4. **Rationale**:
   - These are compounds that BIND to proteins in the complex
   - But are NOT known inhibitors
   - Harder negatives that help the model learn the difference between binding and inhibition

## Implementation in the Improved Pipeline

### Data Files Used

1. **2p2iComplexPairs.txt**
   - Format: `complex_name target_chain target_seq off_target_chain off_target_seq`
   - Contains protein complex information

2. **2p2iInhibitorsSMILES.txt**
   - Format: `complex_name pdb_id complex_id inhibitor_name SMILES label`
   - Contains inhibitor information and SMILES strings

3. **approved_drugs_chemical_structure_identifiers.xlsx**
   - SuperDRUG2 database
   - Contains ~3,633 approved drug compounds

4. **DBD5/** directory
   - Contains PDB files for 282 protein complexes
   - Used for Strategy 2 negative generation

5. **BindersWithComplexname.csv** (pre-processed)
   - Contains binders from BindingDB
   - Already filtered according to the paper's criteria

6. **Binders With Tanimoto Similarity 0.85.csv** (pre-processed)
   - Binders filtered by Tanimoto similarity
   - Only includes compounds dissimilar to known inhibitors

### Preprocessing Steps

#### Step 1: Load Raw Data
```python
complexes_2p2i = load_2p2i_complexes('Data/2p2iComplexPairs.txt')
inhibitors_2p2i = load_2p2i_inhibitors('Data/2p2iInhibitorsSMILES.txt')
```

#### Step 2: Filter Positive Examples
```python
# Keep only complexes with >1 inhibitor
inhibitor_counts = inhibitors_2p2i.groupby('complex_name').size()
complexes_with_multiple_inhibitors = inhibitor_counts[inhibitor_counts > 1].index
positive_examples = inhibitors_2p2i[inhibitors_2p2i['complex_name'].isin(complexes_with_multiple_inhibitors)]
```

#### Step 3: Generate Strategy 1 Negatives
```python
# Combine 2P2I and SuperDRUG2 compounds
all_compounds = list(set(all_2p2i_smiles + superdrug_smiles))

# For each complex, sample random compounds (not known inhibitors)
for complex_name in complexes_list:
    known_inhibitors = get_known_inhibitors(complex_name)
    sample_random_compounds(exclude=known_inhibitors)
```

#### Step 4: Generate Strategy 2 Negatives
```python
# List DBD5 complexes
dbd5_complexes = get_dbd5_complex_list('Data/DBD5/')

# Pair DBD5 complexes with 2P2I inhibitors
for dbd5_complex in dbd5_complexes:
    for inhibitor in random_sample(all_2p2i_inhibitors):
        create_negative_example(dbd5_complex, inhibitor)
```

#### Step 5: Generate Strategy 3 Negatives
```python
# Load pre-processed binders (filtered by BindingDB criteria)
binders_df = pd.read_csv('Data/BindersWithComplexname.csv')
binders_tanimoto = pd.read_csv('Data/Binders With Tanimoto Similarity 0.85.csv')

# Create negative examples from binders
for complex_name in complexes_list:
    complex_binders = get_binders_for_complex(complex_name)
    for binder in complex_binders:
        create_hard_negative(complex_name, binder)
```

#### Step 6: Combine and Save
```python
# Combine all examples
all_examples = pd.concat([
    positive_examples,
    negatives_strategy1,
    negatives_strategy2,
    negatives_strategy3
])

# Save in format compatible with original pipeline
save_to_file('Data/Preprocessed_Dataset_All_Examples.txt')
```

## Output Files

The preprocessing generates these files:

1. **Preprocessed_Dataset_All_Examples.txt**
   - Format: `complex_name target_complex inhibitor_name label`
   - Compatible with the original pipeline
   - Can be used as input to the GNN model

2. **Preprocessed_Dataset_All_Examples.csv**
   - Same data in CSV format
   - Includes additional columns: `strategy`, `smiles`
   - Useful for analysis and debugging

3. **Preprocessed_SMILES_Dict.txt**
   - Maps inhibitor names to SMILES strings
   - Format: `inhibitor_name\tSMILES`

## Dataset Statistics Comparison

| Metric | Research Paper | Expected from Preprocessing |
|--------|---------------|----------------------------|
| **Positive Examples** | 714 | ~700-714 |
| **Negative Strategy 1** | 857 | Variable (10 per complex × ~22 complexes = ~220+) |
| **Negative Strategy 2** | 1,714 | Variable (depends on DBD5 subset used) |
| **Negative Strategy 3** | 11,789 | Depends on pre-processed binders file |
| **Total Negatives** | 10,413 | Variable |
| **Total Examples** | 11,127 | Variable |
| **Unique Complexes** | 22 | ~22 from 2P2I + DBD5 complexes |
| **Pos:Neg Ratio** | 1:14.6 | Similar ratio maintained |

## Key Differences from Original Notebook

### Original `Complete_PPI_Inhibitors_Pipeline_End_To_End.ipynb`:
- Loads pre-processed dataset file: `WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt`
- No preprocessing code included
- Goes directly to model training

### Improved `Improved_PPI_Inhibitors_Pipeline_With_Preprocessing.ipynb`:
- **Includes complete preprocessing pipeline**
- Implements all three negative generation strategies
- Creates dataset from scratch following paper methodology
- Then continues with the same model training as original

## Usage Instructions

### Option 1: Use Pre-processed Dataset (Faster)
If you just want to run the model:
```python
# Use the existing preprocessed file
dataset_file = 'Data/WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt'
```

### Option 2: Regenerate Dataset (Recommended for Research)
If you want to understand or modify the preprocessing:
```python
# Run the improved notebook from the beginning
# It will regenerate the dataset from scratch
# Then use the generated file for training
```

## Validation

To verify your preprocessed dataset matches the paper:

1. **Check positive examples**:
   ```python
   assert len(positive_examples) >= 700, "Should have ~714 positive examples"
   assert positive_examples['complex_name'].nunique() >= 20, "Should have ~22 unique complexes"
   ```

2. **Check negative distribution**:
   ```python
   neg_ratio = len(negatives) / len(positive_examples)
   assert 10 <= neg_ratio <= 20, "Neg:Pos ratio should be ~14.6"
   ```

3. **Check SMILES validity**:
   ```python
   for smiles in all_examples['smiles']:
       mol = Chem.MolFromSmiles(smiles)
       assert mol is not None, f"Invalid SMILES: {smiles}"
   ```

## Citation

If you use this preprocessing pipeline, please cite the original paper:

```bibtex
@article{yaseen2024predicting,
  title={Predicting small-molecule inhibition of protein complexes},
  author={Yaseen, Adiba and Roy, Soumyadip and Akhter, Naeem and Ben-Hur, Asa and Minhas, Fayyaz},
  journal={bioRxiv},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Troubleshooting

### Issue: SuperDRUG2 file not found
**Solution**: Download from the repository or use only 2P2I compounds:
```python
superdrug_smiles = []  # Use empty list if file not available
```

### Issue: Binders file not found
**Solution**: The hard negatives (Strategy 3) require pre-processed binders. You can:
- Skip Strategy 3 (will reduce total negatives)
- Or download binders file from the repository

### Issue: DBD5 PDB files missing
**Solution**: Download DBD5 database or use a subset of available files

## Next Steps

After preprocessing:
1. Load the preprocessed dataset
2. Extract protein features (GNN embeddings, interface features, sequence features)
3. Train the GNN model
4. Perform leave-one-complex-out cross-validation
5. Evaluate on external datasets

See the rest of the `Improved_PPI_Inhibitors_Pipeline_With_Preprocessing.ipynb` for the complete pipeline.
