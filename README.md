# PPI-Inhibitors
2p2iComplexPairs.txt


## Set Up Environment
```
python=3.6
conda install -c conda-forge rdkit
```
## Dataset
2p2iComplexPairs Dataset can be downloaded from this link [https://github.com/adibayaseen/PPI-Inhibitors/blob/01ad4975fb9133825b1bf9e71b64fcdaaa5e4d8b/Data/2p2iComplexPairs.txt]<br/>
2p2iInhibitorsSMILES.txt Dataset can be downloaded from this link [https://github.com/adibayaseen/PPI-Inhibitors/blob/01ad4975fb9133825b1bf9e71b64fcdaaa5e4d8b/Data/2p2iInhibitorsSMILES.txt]<br/>
Binders With Tanimoto Similarity 0.85 dataset from here [Data/Binders With Tanimoto Similarity 0.85.csv] <br/>

# Data discription
## 2p2iComplexPairs.txt
* Input File format <br/>
```
> Name of Complex pair<space> Target chain <space>Target chain Sequence <space>Off Target chain  <space> Off Target chain Sequence<newline> <br/>
1YCR_A_2_B<space>A<space>ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVVvB<space>ETFSDLWKLLPEN <newline><br/>

```
## 2p2iInhibitorsSMILES.txt
* Input File format <br/>
```
> Name of Complex pair<space> Name of the Inhibitor<space>Complex name <space>Inhibted Complex <space> SMILES of the Inhibitor<space>Label<newline> <br/>
1YCR_A_2_B<space>1RV1<space>1YCR<space>IMZ<space>CCOc1cc(ccc1C2=N[C@H]([C@H](N2C(=O)N3CCN(CC3)CCO)c4ccc(cc4)Br)c5ccc(cc5)Br)OC <newline><br/>
1YCR_A_2_B 1RV1 1YCR IMZ CCOc1cc(ccc1C2=N[C@H]([C@H](N2C(=O)N3CCN(CC3)CCO)c4ccc(cc4)Br)c5ccc(cc5)Br)OC  1
1YCR_A_2_B 1T4E 1YCR DIZ c1cc(ccc1[C@H]2C(=O)Nc3ccc(cc3C(=O)N2[C@@H](c4ccc(cc4)Cl)C(=O)O)I)Cl  1
1YCR_A_2_B 2LZG 1YCR 13Q c1cc(cc(c1)Cl)[C@H]2C[C@@H](C(=O)N([C@@H]2c3ccc(cc3)Cl)CC4CC4)CC(=O)O  1
```
## Binders With Tanimoto Similarity 0.85
* Input File format <br/>
```
> (Complex,inhibitor pair)<Next Column> inhibitor Names<Next Column>Binders SMILES<Next Column>Similarity<newline><br/>
1YCR_A_2_B<Next Column> "[array(['IMZ', 'DIZ', '13Q', 'YIN', 'K23', 'MI6', 'TJ2', '07G', 'VZV',
       'LTZ', 'BLF', '0R2', '0R3', '0Y7', 'NUT', '1MN', '1MO', '1MQ',
       '1MT', '1MY', 'Y30', '28W', '2SW', '2TW', '2TZ', '2U0', '2U1',
       '2U5', '2U6', '2U7', '2V8', '35S', '35T', '3UD', '4SS', '4T4',
       '4TH', '62R', '62T', '62Q', 'NUT', '6ZT', '4NJ', '4NX', '6SK',
       '6SJ', '6SS', '6ST', '7HC', '6GG', '6GG', '9QW', 'NUT', 'EYH'],
      dtype='<U3')]"	<Next Column> OC(=O)c1ccc2c3C[C@H]4[C@H]([C@H](c5cccc(Cl)c5F)[C@@]5(N4CC4CC4)C(=O)Nc4nc(Cl)ccc54)n3nc2c1<Next Column> 0.84
 <newline><br/>

```

