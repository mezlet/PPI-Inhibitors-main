"""
PPI Inhibitors - Researcher Dataset Generation Module

This module provides utilities for generating and processing datasets for
protein-protein interaction (PPI) inhibitor prediction. It includes functions
for feature extraction, negative example generation, and dataset loading.

Author: PPI Inhibitors Research Team
"""

import warnings
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, Sampler
from torch.autograd import Variable
import numpy as np
import pandas as pd
import pickle
import math
import random
from itertools import product
from scipy import spatial
from os import listdir
from tqdm import tqdm
import os

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

# BioPython imports
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import PPBuilder, CaPPBuilder

# Sklearn imports
from sklearn.preprocessing import OneHotEncoder, StandardScaler, normalize
from sklearn.model_selection import GroupKFold


# ========================================================================
# CONFIGURATION
# ========================================================================

class Config:
    """Configuration class for paths and parameters"""
    def __init__(self, base_path=None):
        self.base_path = base_path or '/content/PPI-Inhibitors/'
        self.data_path = os.path.join(self.base_path, 'Data/')
        self.features_path = os.path.join(self.base_path, 'Features/')
        self.external_data_path = os.path.join(self.data_path, 'External data/')


# ========================================================================
# DATASET LOADING FUNCTIONS
# ========================================================================

def load_external_dataset(filepath):
    """
    Load external dataset file containing complex IDs, SMILES, and labels.

    Args:
        filepath (str): Path to the dataset file

    Returns:
        list: List of tuples (complex_id, smiles, label)
    """
    with open(filepath) as f:
        lines = f.readlines()

    data = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 3:
            complex_id = parts[0]
            smiles = parts[1]
            label = float(parts[2])
            data.append((complex_id, smiles, label))

    return data


# ========================================================================
# MOLECULAR FINGERPRINT FUNCTIONS
# ========================================================================

def getFP(s, r=3, nBits=2048):
    """
    Generate Morgan fingerprint for a molecule from SMILES string.

    Args:
        s (str): SMILES string
        r (int): Radius for Morgan fingerprint (default: 3)
        nBits (int): Number of bits for fingerprint (default: 2048)

    Returns:
        numpy.ndarray: Fingerprint array or None if invalid SMILES
    """
    compound = Chem.MolFromSmiles(s.strip())
    if compound is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(compound, r, nBits=nBits)
        m = np.zeros((0,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, m)
        return m
    return None


# ========================================================================
# PROTEIN SEQUENCE FEATURE EXTRACTION
# ========================================================================

def twomerFromSeq(s):
    """
    Generate two-mer composition features from protein sequence.
    Groups amino acids into 7 categories and counts all 2-mer combinations.

    Args:
        s (str): Protein sequence string

    Returns:
        numpy.ndarray: Array of 49 features (7x7 combinations)
    """
    k = 2
    groups = {
        'A':'1', 'V':'1', 'G':'1',
        'I':'2', 'L':'2', 'F':'2', 'P':'2',
        'Y':'3', 'M':'3', 'T':'3', 'S':'3',
        'H':'4', 'N':'4', 'Q':'4', 'W':'4',
        'R':'5', 'K':'5',
        'D':'6', 'E':'6',
        'C':'7'
    }

    crossproduct = [''.join(i) for i in product("1234567", repeat=k)]
    for i in range(0, len(crossproduct)):
        crossproduct[i] = int(crossproduct[i])

    ind = []
    for i in range(0, len(crossproduct)):
        ind.append(i)

    combinations = dict(zip(crossproduct, ind))
    V = np.zeros(int((math.pow(7, k))))  # vector of 49 length with zero entries

    try:
        for j in range(0, len(s) - k + 1):
            kmer = s[j:j+k]
            c = ''
            for l in range(0, k):
                c += groups[kmer[l]]
            V[combinations[int(c)]] += 1
    except:
        count = {'1':0, '2':0, '3':0, '4':0, '5':0, '6':0, '7':0}
        for q in range(0, len(s)):
            if s[q] == 'A' or s[q] == 'V' or s[q] == 'G':
                count['1'] += 1
            if s[q] == 'I' or s[q] == 'L' or s[q] == 'F' or s[q] == 'P':
                count['2'] += 1
            if s[q] == 'Y' or s[q] == 'M' or s[q] == 'T' or s[q] == 'S':
                count['3'] += 1
            if s[q] == 'H' or s[q] == 'N' or s[q] == 'Q' or s[q] == 'W':
                count['4'] += 1
            if s[q] == 'R' or s[q] == 'K':
                count['5'] += 1
            if s[q] == 'D' or s[q] == 'E':
                count['6'] += 1
            if s[q] == 'C':
                count['7'] += 1

        val = list(count.values())
        key = list(count.keys())
        m = 0
        ind = 0
        for t in range(0, len(val)):
            if m < val[t]:
                m = val[t]
                ind = t
        m = key[ind]  # group number of maximum occurring group alphabets

        for j in range(0, len(s) - k + 1):
            kmer = s[j:j+k]
            c = ''
            for l in range(0, k):
                if kmer[l] not in groups:
                    c += m
                else:
                    c += groups[kmer[l]]
            V[combinations[int(c)]] += 1

    V = V / (len(s) - 1)
    return np.array(V)


def prot_feats_seq(seq):
    """
    Extract protein sequence features including:
    - Amino acid composition (20 features, normalized)
    - Two-mer composition (49 features, normalized)

    Args:
        seq (str): Protein sequence string

    Returns:
        numpy.ndarray: Array of 69 features (20 + 49)
    """
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    f = []

    X = ProteinAnalysis(str(seq))
    X.molecular_weight()  # throws an error if 'X' in sequence. we skip such sequences
    p = X.get_amino_acids_percent()

    dp = []
    for a in aa:
        dp.append(p[a])
    dp = np.array(dp)
    dp = normalize(np.atleast_2d(dp), norm='l2', copy=True, axis=1, return_norm=False)
    f.extend(dp[0])

    tm = np.array(twomerFromSeq(str(seq)))
    tm = normalize(np.atleast_2d(tm), norm='l2', copy=True, axis=1, return_norm=False)
    f.extend(tm[0])

    return np.array(f)


# ========================================================================
# PROTEIN-PROTEIN INTERFACE FEATURES
# ========================================================================

def make_dic():
    """
    Create a dictionary with all possible amino acid pairs as keys,
    initialized to 0.0.

    Returns:
        dict: Dictionary with amino acid pair tuples as keys
    """
    prot_dic = {}
    letters = IUPACData.protein_letters
    for i in range(len(letters)):
        for j in range(i, len(letters)):
            prot_dic[(letters[i], letters[j])] = 0.0
    prot_dic[('_', '_')] = 0.0  # for amino acids other than 20 natural
    return prot_dic


def generate_pair_features(dist_info, xl, xr):
    """
    Generate amino acid pair interaction features from distance information.
    Counts occurrences of each amino acid pair within interaction distance.

    Args:
        dist_info (list): List of tuples with distance information
        xl (list): Left chain residues
        xr (list): Right chain residues

    Returns:
        dict: Dictionary with amino acid pair counts
    """
    prot_dic = make_dic()

    for rec in dist_info:
        try:
            l_letter = three_to_one(xl[rec[0]].get_resname())
            r_letter = three_to_one(xr[rec[1]].get_resname())

            if (l_letter, r_letter) in prot_dic.keys():
                prot_dic[(l_letter, r_letter)] += 1
            elif (r_letter, l_letter) in prot_dic.keys():
                prot_dic[(r_letter, l_letter)] += 1
        except:
            prot_dic[('_', '_')] += 1

    return prot_dic


def getCoords(R):
    """
    Get atom coordinates given a list of biopython residues.

    Args:
        R (list): List of BioPython residue objects

    Returns:
        list: List of coordinate arrays for each residue
    """
    Coords = []
    for (idx, r) in enumerate(R):
        v = [ak.get_coord() for ak in r.get_list()]
        Coords.append(v)
    return Coords


def getDist(C0, C1, thr=np.inf):
    """
    Calculate pairwise distances between two chains.

    Args:
        C0 (list): Coordinates of first chain
        C1 (list): Coordinates of second chain
        thr (float): Distance threshold (default: infinity)

    Returns:
        tuple: (N0, N1) with neighbor lists
    """
    N0 = []
    N1 = []
    for i in range(len(C0)):
        for j in range(len(C1)):
            d = spatial.distance.cdist(C0[i], C1[j]).min()
            if (d < thr):
                N0.append((i, j, d))
                N1.append((j, i, d))
    return (N0, N1)


def extract_feats(dic, key_list_path=None):
    """
    Extract features from amino acid pair dictionary in consistent order.
    Uses pre-saved key list to ensure features are in the same order.

    Args:
        dic (dict): Dictionary with amino acid pair counts
        key_list_path (str): Path to key list file (optional)

    Returns:
        list: List of feature values in consistent order
    """
    feats = []
    if key_list_path is None:
        key_list_path = '/content/PPI-Inhibitors/Features/prote_letter_pair_keys.npy'

    key_list = np.load(key_list_path)
    for key in key_list:
        feats.append(dic[(key[0].decode('utf-8'), key[1].decode('utf-8'))])
    return feats


def chainLabel(Cname_T, xl_T, Cname, xl):
    """
    Generate interface features between two protein chains.

    Args:
        Cname_T (str): Target chain name
        xl_T (list): Target chain residues
        Cname (str): Off-target chain name
        xl (list): Off-target chain residues

    Returns:
        numpy.ndarray: Array of interface features (amino acid pair counts)
    """
    tc = getCoords(xl_T)
    nc = getCoords(xl)
    D = getDist(tc, nc, thr=8.0)
    feats = extract_feats(generate_pair_features(D, xl_T, xl))
    return feats


def Struct2chain(stx):
    """
    Extract chains from PDB structure with sequence information.

    Args:
        stx (str): Path to PDB file

    Returns:
        list: List of tuples (chain_name, sequence, sequence_length, residue_list)
    """
    p = PDBParser()
    L = []
    stx = p.get_structure('X', stx)

    for model in stx:
        for C in model:
            RL = []
            for R in C:
                RL.append(R)

            pp = PPBuilder().build_peptides(C)
            if len(pp) == 0:
                pp = CaPPBuilder().build_peptides(C)

            seq = ''.join([str(p.get_sequence()) for p in pp])
            seq_L = len(seq)
            L.append((C.full_id[2], seq, seq_L, RL))

    return L


def InterfaceFeatures(Complexs, pdbloc):
    """
    Calculate interface features for protein complexes.
    Computes amino acid pair interaction features at protein-protein interfaces.

    Args:
        Complexs (list): List of complex names
        pdbloc (str): Path to directory containing PDB files

    Returns:
        dict: Dictionary mapping complex names to interface features
    """
    Found = listdir(pdbloc)
    InterfaceFeatures_dict = []
    InterfaceFeatures_dict = dict(InterfaceFeatures_dict)
    comp_id = list(set(Complexs))

    for ids in range(len(comp_id)):
        if comp_id[ids] + '.pdb' in Found:
            stx = pdbloc + '/' + comp_id[ids] + '.pdb'
            chains = Struct2chain(stx)

            for j in range(len(chains)):
                Cname_T, seq_T, L_T, xl_T = chains[j]
                for k in range(j, len(chains)):
                    Cname, seq, L, xl = chains[k]
                    name = comp_id[ids]
                    Interface = chainLabel(Cname_T, xl_T, Cname, xl)
                    InterfaceF = np.array(Interface)
                    InterfaceF = normalize(np.atleast_2d(InterfaceF), norm='l2',
                                         copy=True, axis=1, return_norm=False)
                    if name not in InterfaceFeatures_dict.keys():
                        InterfaceFeatures_dict[name] = Interface

    return InterfaceFeatures_dict


# ========================================================================
# GRAPH NEURAL NETWORK DATA PROCESSING
# ========================================================================

def atom1(structure):
    """
    One-hot encode atom types (13 types).

    Args:
        structure: BioPython structure object

    Returns:
        numpy.ndarray: One-hot encoded atom types
    """
    atomslist = np.array(sorted(['C', 'CA', 'CB', 'CG', 'CH2', 'N', 'NH2',
                                  'OG', 'OH', 'O1', 'O2', 'SE', '1'])).reshape(-1, 1)
    enc = OneHotEncoder(handle_unknown='ignore')
    enc.fit(atomslist)

    atom_list = []
    for atom in structure.get_atoms():
        if atom.get_name() in atomslist:
            atom_list.append(atom.get_name())
        else:
            atom_list.append("1")

    atoms_onehot = enc.transform(np.array(atom_list).reshape(-1, 1)).toarray()
    return atoms_onehot


def res1(structure):
    """
    One-hot encode residue types (21 types).

    Args:
        structure: BioPython structure object

    Returns:
        numpy.ndarray: One-hot encoded residue types
    """
    residuelist = np.array(sorted(['ALA', 'ARG', 'ASN', 'ASP', 'GLN', 'GLU',
                                    'GLY', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                                    'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                                    'CYS', 'HIS', '1'])).reshape(-1, 1)
    encr = OneHotEncoder(handle_unknown='ignore')
    encr.fit(residuelist)

    residue_list = []
    for atom in structure.get_atoms():
        if atom.get_parent().get_resname() in residuelist:
            residue_list.append(atom.get_parent().get_resname())
        else:
            residue_list.append("1")

    res_onehot = encr.transform(np.array(residue_list).reshape(-1, 1)).toarray()
    return res_onehot


def neigh1(structure, cutoff=6.0, max_neighbors=10):
    """
    Calculate neighbors for each atom (same residue and different residue).

    Args:
        structure: BioPython structure object
        cutoff (float): Distance cutoff for neighbors (default: 6.0 Angstroms)
        max_neighbors (int): Maximum number of neighbors to store (default: 10)

    Returns:
        tuple: (neigh_same_res, neigh_diff_res) neighbor arrays
    """
    atom_list = np.array([atom for atom in structure.get_atoms()])
    p4 = NeighborSearch(atom_list)
    neighbour_list = p4.search_all(cutoff, level="A")
    neighbour_list = np.array(neighbour_list)

    dist = np.array([nl[0] - nl[1] for nl in neighbour_list])
    place = np.argsort(dist)
    sorted_neighbour_list = neighbour_list[place]

    source_vertex_list = np.array(sorted_neighbour_list[:, 0])
    neighbour_vertex_list = np.array(sorted_neighbour_list[:, 1])

    old_atom_number = [atom.get_serial_number() for atom in atom_list]
    old_residue_number = [atom.get_parent().get_id()[1] for atom in atom_list]
    old_atom_number = np.array(old_atom_number)
    old_residue_number = np.array(old_residue_number)

    total_atoms = len(atom_list)
    neigh_same_res = np.array([[-1] * max_neighbors for _ in range(total_atoms)])
    neigh_diff_res = np.array([[-1] * max_neighbors for _ in range(total_atoms)])
    same_flag = [0] * total_atoms
    diff_flag = [0] * total_atoms

    for i in range(len(source_vertex_list)):
        source_atom_id = source_vertex_list[i].get_serial_number()
        neigh_atom_id = neighbour_vertex_list[i].get_serial_number()
        source_atom_res = source_vertex_list[i].get_parent().get_id()[1]
        neigh_atom_res = neighbour_vertex_list[i].get_parent().get_id()[1]

        temp_index1 = np.where(source_atom_id == old_atom_number)[0]
        temp_index2 = np.where(neigh_atom_id == old_atom_number)[0]

        source_index = None
        neigh_index = None

        for i1 in temp_index1:
            if old_residue_number[i1] == source_atom_res:
                source_index = i1
                break

        for i1 in temp_index2:
            if old_residue_number[i1] == neigh_atom_res:
                neigh_index = i1
                break

        if source_index is None or neigh_index is None:
            continue

        if source_atom_res == neigh_atom_res:
            if same_flag[source_index] < max_neighbors:
                neigh_same_res[source_index][same_flag[source_index]] = neigh_index
                same_flag[source_index] += 1
            if same_flag[neigh_index] < max_neighbors:
                neigh_same_res[neigh_index][same_flag[neigh_index]] = source_index
                same_flag[neigh_index] += 1
        else:
            if diff_flag[source_index] < max_neighbors:
                neigh_diff_res[source_index][diff_flag[source_index]] = neigh_index
                diff_flag[source_index] += 1
            if diff_flag[neigh_index] < max_neighbors:
                neigh_diff_res[neigh_index][diff_flag[neigh_index]] = source_index
                diff_flag[neigh_index] += 1

    return neigh_same_res, neigh_diff_res


def processProtein(UniqueProtein, PdBloc, device=None):
    """
    Process PDB files to create graph neural network data.
    Converts protein structures into graph representations with:
    - Node features: atom types and residue types (one-hot encoded)
    - Edge features: neighbor lists (same/different residue)

    Args:
        UniqueProtein (list): List of unique protein PDB names
        PdBloc (str): Path to directory containing PDB files
        device (torch.device): Device to use (CPU/GPU)

    Returns:
        dict: Dictionary mapping protein names to GNN data
    """
    if device is None:
        use_cuda = torch.cuda.is_available()
        device = torch.device("cuda:0" if use_cuda else "cpu")

    PData_dict = {}

    for i in range(len(UniqueProtein)):
        UniqueProtein[i] = UniqueProtein[i].split('.pdb')[0]
        P1 = PdBloc + UniqueProtein[i] + '.pdb'

        parser = PDBParser()
        with warnings.catch_warnings(record=True) as w:
            structure = parser.get_structure("", P1)

        one_hot_atom = (atom1(structure))
        one_hot_res = (res1(structure))
        neigh_same_res, neigh_diff_res = (neigh1(structure))

        # Make the graph ready for PyTorch Geometric GCN algorithms:
        one_hot_atom = torch.tensor(one_hot_atom, dtype=torch.float32).to(device)
        one_hot_res = torch.tensor(one_hot_res, dtype=torch.float32).to(device)
        neigh_same_res = torch.tensor(neigh_same_res).to(device).long()
        neigh_diff_res = torch.tensor(neigh_diff_res).to(device).long()

        GNNData = [one_hot_atom, one_hot_res, neigh_same_res, neigh_diff_res]
        PData_dict[UniqueProtein[i]] = GNNData

    return PData_dict


# ========================================================================
# LOAD PROTEIN FEATURES (SVM + INTERFACE)
# ========================================================================

def LoadProtein_SVM_Features(UniqueProtein, Pdbloc):
    """
    Load and compute SVM features for proteins.
    Combines sequence features and interface features.

    Args:
        UniqueProtein (list): List of unique protein names
        Pdbloc (str): Path to directory containing PDB files

    Returns:
        tuple: (InterfaceFeatures, SequenceFeatures, AllFeatures) dictionaries
    """
    pdbname = listdir(Pdbloc)
    InterfaceFeatures_dict = []
    InterfaceFeatures_dict = dict(InterfaceFeatures_dict)
    SequenceFeatures_dict = []
    SequenceFeatures_dict = dict(SequenceFeatures_dict)
    AllFeatures_dict = []
    AllFeatures_dict = dict(AllFeatures_dict)

    for b in range(len(UniqueProtein)):
        if UniqueProtein[b] + '.pdb' in pdbname:
            stx = Pdbloc + UniqueProtein[b] + '.pdb'
            chains = Struct2chain(stx)

            # Interface Features
            for j in range(len(chains)):
                Cname_T, seq_T, L_T, xl_T = chains[j]
                for k in range(j, len(chains)):
                    Cname, seq, L, xl = chains[k]
                    name = UniqueProtein[b]
                    Interface = chainLabel(Cname_T, xl_T, Cname, xl)
                    seq_TF = prot_feats_seq(seq_T)
                    seq_NTF = prot_feats_seq(seq)
                    SeQFeatures = (seq_TF + seq_NTF) / 2
                    InterfaceF = np.array(Interface)
                    InterfaceF = normalize(np.atleast_2d(InterfaceF), norm='l2',
                                         copy=True, axis=1, return_norm=False)
                    if name not in InterfaceFeatures_dict.keys():
                        InterfaceFeatures_dict[name] = Interface
                        SequenceFeatures_dict[name] = SeQFeatures
                        AllFeatures_dict[name] = np.append(SeQFeatures, Interface)

    return InterfaceFeatures_dict, SequenceFeatures_dict, AllFeatures_dict


# ========================================================================
# NEGATIVE EXAMPLE GENERATION
# ========================================================================

def External_GenerateRandomNegative(posexamples, config=None):
    """
    Generate random negative examples for training.
    Uses two strategies:
    1. Pair known complexes with random drugs from SuperDrug database
    2. Pair known ligands with random complexes from Ubench5

    Args:
        posexamples (dict): Dictionary of positive examples (complex, ligand) -> (complex, smiles)
        config (Config): Configuration object with paths

    Returns:
        tuple: (AllPos, AllNeg, SuperDrug_dict)
    """
    if config is None:
        config = Config()

    NegtiveRatio = 1

    # SuperDrugbank Names
    superdrug_path = os.path.join(config.data_path,
                                   'approved_drugs_chemical_structure_identifiers.xlsx')
    SuperdrugNames = pd.read_excel(superdrug_path, usecols="B").values
    SuperdrugNames = SuperdrugNames[1:]
    SuperdrugNames = np.array([s[0] for s in SuperdrugNames])

    # SuperDrugbank SMILES
    df_Superdrug = pd.read_excel(superdrug_path, usecols="C").values
    df_Superdrug = df_Superdrug[1:]
    df_Superdrug_Compounds = np.array([c[0] for c in df_Superdrug])
    SuperDrug_dict = dict(zip(SuperdrugNames, df_Superdrug_Compounds))

    # Load DBD5 protein data
    dbd5_path = os.path.join(config.features_path, 'NewUbench5InterfaceandSeq_dict.npy')
    DBD5_ProteinData_dict = pickle.load(open(dbd5_path, "rb"))
    Ubench5CompNames = list(set(list(DBD5_ProteinData_dict.keys())))

    AllNeg = []
    AllPos = []
    complex_ligand_dict = {}

    for key, val in posexamples:
        if key not in complex_ligand_dict:
            complex_ligand_dict[key] = posexamples[key, val][1]
        else:
            complex_ligand_dict[key] = np.append(
                complex_ligand_dict.get(key, ()),
                posexamples[key, val][1]
            )

    Complexnames = list(complex_ligand_dict.keys())
    totalcomp = list(set(complex_ligand_dict.keys()))

    # Method 1: Pair complexes with random drugs
    for everycomp in totalcomp:
        origanlL = complex_ligand_dict[everycomp]
        pos = [(everycomp, origanlL[t]) for t in range(len(origanlL))]
        NN = NegtiveRatio * len(pos)
        negs = []
        AllPos.extend(pos)

        while (len(negs) < NN):
            LigandR = random.choice(SuperdrugNames)
            LigandR_smile = SuperDrug_dict[LigandR]
            Npair = ((everycomp, LigandR_smile))
            if (LigandR not in origanlL and
                Npair not in AllNeg and
                Npair not in AllPos and
                getFP(LigandR_smile) is not None):
                negs.append(Npair)

        AllNeg.extend(negs)

    # Method 2: Cross-reference with other complexes
    for everycomp in totalcomp:
        origanlL = complex_ligand_dict[everycomp]
        pos = [(everycomp, origanlL[t]) for t in range(len(origanlL))]
        NN = NegtiveRatio * len(pos)
        negs = []

        while (len(negs) < NN):
            for everyL in origanlL:
                ComplexR = random.choice(Ubench5CompNames)
                Npair = ((ComplexR, everyL))
                if (ComplexR != everycomp and
                    Npair not in AllNeg and
                    Npair not in AllPos):
                    negs.append(Npair)

        AllNeg.extend(negs)

    return np.array(AllPos), np.array(AllNeg), SuperDrug_dict


# ========================================================================
# PREDICTION FUNCTIONS
# ========================================================================

def PredictScorefromFile(filename, Pdbloc, Pscaler, Cscaler, trainedModel_IPPI,
                         train_GNN, LOCOcomplexname, config=None, device=None):
    """
    Predict interaction scores for protein-ligand pairs from file.

    Args:
        filename (str): Path to input file with format: PdbId InhibitedComplex LigandId SMILES
        Pdbloc (str): Path to PDB files directory
        Pscaler (StandardScaler): StandardScaler for protein features
        Cscaler (StandardScaler): StandardScaler for compound features
        trainedModel_IPPI: Trained MLP model
        train_GNN: Trained GNN model
        LOCOcomplexname (str): Leave-one-complex-out name
        config (Config): Configuration object
        device (torch.device): Device to use (CPU/GPU)

    Returns:
        tuple: (predictions, targets)
    """
    if config is None:
        config = Config()

    if device is None:
        use_cuda = torch.cuda.is_available()
        device = torch.device("cuda:0" if use_cuda else "cpu")

    with open(filename) as f:
        D = f.readlines()

    InhibitedComp = []
    PdbId = []
    Ligandnames = []
    SMILES = []
    targets = []
    All_data_list = []

    # Parse input file
    for d in tqdm(D):
        Pdbid, inhibtedc, Ligandid, smiles = d.split()
        if getFP(smiles) is not None:
            PdbId.append(Pdbid)
            Ligandnames.append(Ligandid)
            SMILES.append(smiles)
            InhibitedComp.append(inhibtedc)

    # Load pre-generated examples with labels
    complexnames = []
    SMILES = []
    targets = []

    external_file = os.path.join(config.external_data_path,
                                  filename.split('.txt')[0].split('/')[-1] +
                                  '_External_All_Examples.txt')
    with open(external_file) as f:
        D = f.readlines()

    for d in tqdm(D):
        complexname, smiles, target = d.split()
        complexnames.append(complexname)
        SMILES.append(smiles)
        targets.append(target)

    # Get unique proteins and process them
    pdbname = listdir(Pdbloc)
    mypdb = []
    for p in pdbname:
        if p.split('.pdb')[0] in PdbId:
            mypdb.append(p)

    UniqueProtein = list(set(mypdb))
    External_Protein_GNN_Data_dict = processProtein(UniqueProtein, Pdbloc, device)

    # Get Sequence + Interface features
    s, i, External_ProteinSeqandInterfaceData_dict = LoadProtein_SVM_Features(
        UniqueProtein, Pdbloc
    )

    # Merge with DBD5 features
    ubench5_path = os.path.join(config.features_path, 'NewUbench5InterfaceandSeq_dict.npy')
    Ubench5InterfaceandSeq_dict = pickle.load(open(ubench5_path, "rb"))
    All_External_ProteinSeqandInterfaceData_dict = dict(
        list(External_ProteinSeqandInterfaceData_dict.items()) +
        list(Ubench5InterfaceandSeq_dict.items())
    )

    # Load DBD5 GNN data
    dbd5_gnn_path = os.path.join(config.features_path, 'DBD5_ProteinData_dict.pickle')
    DBD5_Protein_GNN_Data_dict = pickle.load(open(dbd5_gnn_path, "rb"))
    All_Protein_GNN_Data_dict = dict(
        list(External_Protein_GNN_Data_dict.items()) +
        list(DBD5_Protein_GNN_Data_dict.items())
    )

    # Move GNN data to device
    for d in All_Protein_GNN_Data_dict:
        data = All_Protein_GNN_Data_dict[d]
        All_Protein_GNN_Data_dict[d] = [
            data[0].to(device), data[1].to(device),
            data[2].to(device), data[3].to(device)
        ]

    # Prepare test data
    Cttname = []
    Ctt = []
    Pttname = []
    Ptt = []

    for (complexname, ligandsmile) in zip(complexnames, SMILES):
        Cttname.append(ligandsmile)
        Ctt.append(getFP(ligandsmile))
        Pttname.append(complexname)
        Ptt.append(All_External_ProteinSeqandInterfaceData_dict[complexname])

    # Standardization
    Ctt = Cscaler.transform(Ctt)
    Cttdict = dict(zip(Cttname, torch.FloatTensor(Ctt).to(device)))
    Ptt = Pscaler.transform(Ptt)
    Pttdict = dict(zip(Pttname, torch.FloatTensor(Ptt).to(device)))

    # Prediction
    Y_t, Z, Targets = [], [], []

    # Set models to evaluation mode
    trainedModel_IPPI.eval()
    train_GNN.eval()

    with torch.no_grad():
        for target, (complexname, ligandsmile) in zip(targets, zip(complexnames, SMILES)):
            # 1. Get GNN features (already 2D: [1, 512])
            GNN_features = train_GNN(All_Protein_GNN_Data_dict[complexname])

            # 2. Get Ligand and Interface features (which are 1D)
            ligand_feats_1d = Cttdict[ligandsmile]
            interface_feats_1d = Pttdict[complexname]

            # 3. Add a batch dimension (unsqueeze) to make them 2D
            ligand_feats_2d = ligand_feats_1d.unsqueeze(0)       # Shape becomes [1, 2048]
            interface_feats_2d = interface_feats_1d.unsqueeze(0)  # Shape becomes [1, 280]

            # 4. Pass all 2D tensors to the MLP
            test_score = trainedModel_IPPI(GNN_features, ligand_feats_2d, interface_feats_2d)

            test_score = test_score.cpu().data.numpy()[0]
            Z.append(test_score)
            Targets.append(float(target))

    return Z, Targets


# ========================================================================
# UTILITY FUNCTIONS
# ========================================================================

def get_device():
    """
    Get the best available device (CUDA or CPU).

    Returns:
        torch.device: Device to use for computations
    """
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")
    return device


if __name__ == "__main__":
    print("PPI Inhibitors - Researcher Dataset Generation Module")
    print("=" * 60)
    print("This module provides utilities for generating and processing")
    print("datasets for protein-protein interaction inhibitor prediction.")
    print("\nAvailable functions:")
    print("  - load_external_dataset()")
    print("  - getFP()")
    print("  - processProtein()")
    print("  - LoadProtein_SVM_Features()")
    print("  - External_GenerateRandomNegative()")
    print("  - PredictScorefromFile()")
    print("\nFor detailed documentation, see function docstrings.")
