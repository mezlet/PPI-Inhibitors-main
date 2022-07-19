# PPI-Inhibitors
2p2iComplexPairs.txt
2p2iInhibitorsSMILES.txt
Binders With Tanimoto Similarity 0.85

## Set Up Environment
```
python=3.6
conda install -c conda-forge rdkit
```
## Dataset
2p2iComplexPairs Dataset can be downloaded from this link [https://github.com/adibayaseen/HKRCPI/tree/main/Datasets/NR-HCPI]<br/>
2p2iInhibitorsSMILES.txt Dataset can be downloaded from this link [https://github.com/adibayaseen/HKRCPI/tree/main/Datasets/NR-HCPI]<br/>
Binders With Tanimoto Similarity 0.85 dataset from here https://github.com/adibayaseen/HKRCPI/blob/main/Datasets/BindingDB_m62021_top4000_1783000nM.txt <br/>
SuperDrugbank2 dataset from here https://github.com/adibayaseen/HKRCPI/blob/main/Datasets/approved_drugs_chemical_structure_identifiers.xlsx <br/>

# Data discription
## 2p2iComplexPairs.txt
* Input File format <br/>
```
> Name of Complex pair<space> Target chain <space>Target chain Sequence <space>Off Target chain  <space> Off Target chain Sequence<newline> <br/>
1YCR_A_2_B<space>A<space>ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVVvB<space>ETFSDLWKLLPEN <newline><br/>

```

