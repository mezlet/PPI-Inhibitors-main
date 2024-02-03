# Predicting small-molecule inhibition of protein complexes

## Abstract

**Motivation:**
<br>Protein-Protein Interactions (PPIs) have key roles in biological processes and are involved in numerous diseases, making discovery of PPI inhibitors a key focus in drug development. Traditional approaches for this purpose involve expensive and time-consuming experiments in the wet lab. Such drug development efforts can greatly benefit from novel machine learning methods to predict inhibitors of specific PPIs. To the best of our knowledge, there is no existing method that takes a protein complex and a compound as inputs to predict whether the given compound acts as an inhibitor of that complex.</br>
**Methods:**
<br>We present a graph neural network approach that takes the structure of a protein complex and the SMILES representation of a compound to predict the potential of the compound to inhibit the interaction between proteins in the given complex in a targeted manner.</br>
**Results:**
<br>The proposed method, trained and validated on 714 inhibitors of 23 complexes in the 2p2i-DB-v2 database, demonstrates superior predictive accuracy (cross-validation AUC-ROC of 0.85) outperforming baseline kernel methods and pre-trained neural network approaches. We further tested the predictive performance of our model on two independent external datasets – one collected from recent publications and another consisting of putative inhibitors of the SARS-CoV-2-Spike and Human-ACE2 protein complex with AUC-ROCs of 0.82 and 0.78, respectively. This pioneering application represents the first successful method for predicting specific protein complex inhibitors for drug design and development and lays the groundwork for future model development in this vital field.</br>

## Set Up Environment
```
python=3.6
conda install -c conda-forge rdkit
```
## Dataset
2p2iComplexPairs Dataset can be downloaded from this link [https://github.com/adibayaseen/PPI-Inhibitors/blob/01ad4975fb9133825b1bf9e71b64fcdaaa5e4d8b/Data/2p2iComplexPairs.txt]<br/>
2p2iInhibitorsSMILES.txt Dataset can be downloaded from this link [https://github.com/adibayaseen/PPI-Inhibitors/blob/01ad4975fb9133825b1bf9e71b64fcdaaa5e4d8b/Data/2p2iInhibitorsSMILES.txt]<br/>
Binders With Tanimoto Similarity 0.85 dataset from here [https://github.com/adibayaseen/PPI-Inhibitors/blob/01ad4975fb9133825b1bf9e71b64fcdaaa5e4d8b/Data/Binders%20With%20Tanimoto%20Similarity%200.85.csv] <br/>

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
GNN-based pipeline Protein complex Features for Positive data (2p2i complexes) Drive link is here <br> [https://drive.google.com/file/d/1goeDiPZSKT1Xx3j00eNG9xlqYkLLv1gW/view?usp=sharing]<br/>
GNN-based pipeline DBD5 Protein complex Features for Negative examples, Drive link is here  <br> [https://drive.google.com/file/d/1GOYEKLQCoGea9QQ72kujy0rdJKbUSYAE/view?usp=sharing]<br/>
# Code
## Setting the path of Input 
 <br> Clone GitHub repository to collab and give githubpath as /content/PPI-Inhibitors </br>
 <br> Download 2p2i complex features and DBD5 complex features from Google Drive and path with folder name GNN-PPI-Inhibitor  </br>
## GNN-based Results
GNN-based pipeline Result can be reconstructed using this collab code <br> [https://github.com/adibayaseen/PPI-Inhibitors/blob/a0c8476b075d6da82bcde804db09502dfd8d0100/code/GNN_based_Pipeline_Load_Model_ipynb.ipynb] </br>
##  GearNet Embedding in combination with GNN-based pipeline Results
GearNet Embedding in combination with GNN-based pipeline can be reconstructed using the link <br> [https://github.com/adibayaseen/PPI-Inhibitors/blob/58a1adced4c0efede92d5689ecf22ad2755bd03f/code/Load_Model_GearNetEmbedding_and_GNN_based_pipeline_ipynb.ipynb]</br>

