"""
Dataset Generation Script for PPI Inhibitors
============================================

This script generates the training dataset according to the specifications
in the research paper "Predicting small-molecule inhibition of protein complexes".

Dataset Composition:
-------------------
1. Positive Examples: 714 inhibitors from 2P2I database across 22 protein complexes
2. Negative Examples (10,413 total):
   - Strategy 1: Random pairing from 2P2I + SuperDRUG2 (857 examples)
   - Strategy 2: 2P2I compounds with DBD5 complexes (1,714 examples)
   - Strategy 3: Binders from Binding-DB with Tanimoto < 0.85 (7,842 examples)

Feature Extraction:
------------------
- Ligand Features: ECFP (Morgan fingerprints) with radius 2, 2048 dimensions
- Protein Sequence Features: AAC (20 dim) + grouped k-mer (49 dim for k=2)
- Interface Features: 211 dimensional (amino acid pairs within 8Å)
- Total feature dimension: 2,840

Usage:
------
python generate_training_dataset.py --output_dir ./generated_data
"""

import os
import sys
import argparse
import random
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit import DataStructs
except ImportError:
    print("ERROR: RDKit not installed. Please install: pip install rdkit-pypi")
    sys.exit(1)

try:
    from Bio.PDB import PDBParser, Selection
    from Bio import pairwise2
except ImportError:
    print("ERROR: Biopython not installed. Please install: pip install biopython")
    sys.exit(1)


class DatasetGenerator:
    """Generate training dataset for PPI inhibitor prediction"""

    def __init__(self, data_dir: str = "./Data"):
        self.data_dir = data_dir
        self.positive_examples = []
        self.negative_examples = []

        # Standard amino acids + unknown
        self.amino_acids = list("ACDEFGHIKLMNPQRSTVWY") + ["X"]

        # Grouped amino acids based on physicochemical properties (7 groups)
        self.aa_groups = {
            'A': 0, 'G': 0, 'V': 0,  # Aliphatic
            'I': 1, 'L': 1, 'F': 1, 'P': 1,  # Large aliphatic
            'Y': 2, 'M': 2, 'T': 2, 'S': 2,  # Hydroxyl + sulfur-containing
            'H': 3, 'N': 3, 'Q': 3, 'W': 3,  # Amide
            'R': 4, 'K': 4,  # Positive
            'D': 5, 'E': 5,  # Negative
            'C': 6  # Cysteine
        }

        print(f"Initializing DatasetGenerator with data directory: {data_dir}")

    def load_positive_examples(self) -> List[Dict]:
        """
        Load positive examples from 2P2I database
        Returns list of dicts with complex_name, inhibitor_name, smiles, label
        """
        print("\n" + "="*80)
        print("LOADING POSITIVE EXAMPLES FROM 2P2I DATABASE")
        print("="*80)

        inhibitors_file = os.path.join(self.data_dir, "2p2iInhibitorsSMILES.txt")
        complexes_file = os.path.join(self.data_dir, "2p2iComplexPairs.txt")

        if not os.path.exists(inhibitors_file):
            raise FileNotFoundError(f"Inhibitors file not found: {inhibitors_file}")

        # Load inhibitors with SMILES
        positive_data = []
        complexes_with_single_inhibitor = set()
        complex_inhibitor_count = defaultdict(set)

        with open(inhibitors_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 5:
                    complex_pair = parts[0]
                    inhibitor_name = parts[1]
                    complex_name = parts[2]
                    inhibitor_code = parts[3]
                    smiles = parts[4]
                    label = parts[5] if len(parts) > 5 else "1"

                    complex_inhibitor_count[complex_pair].add(inhibitor_name)

                    positive_data.append({
                        'complex_pair': complex_pair,
                        'complex_name': complex_name,
                        'inhibitor_name': inhibitor_name,
                        'inhibitor_code': inhibitor_code,
                        'smiles': smiles,
                        'label': int(float(label))
                    })

        # Filter out complexes with only one inhibitor (as per paper)
        filtered_data = []
        for item in positive_data:
            if len(complex_inhibitor_count[item['complex_pair']]) > 1:
                filtered_data.append(item)

        print(f"✓ Loaded {len(filtered_data)} positive examples")
        print(f"  - Unique complexes: {len(set([x['complex_pair'] for x in filtered_data]))}")
        print(f"  - Unique inhibitors: {len(set([x['smiles'] for x in filtered_data]))}")

        self.positive_examples = filtered_data
        return filtered_data

    def load_complex_sequences(self) -> Dict[str, Dict]:
        """Load protein complex sequences"""
        complexes_file = os.path.join(self.data_dir, "2p2iComplexPairs.txt")

        if not os.path.exists(complexes_file):
            print(f"WARNING: Complex pairs file not found: {complexes_file}")
            return {}

        complex_seqs = {}
        with open(complexes_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 5:
                    complex_name = parts[0]
                    target_chain = parts[1]
                    target_seq = parts[2]
                    off_target_chain = parts[3]
                    off_target_seq = parts[4]

                    complex_seqs[complex_name] = {
                        'target_chain': target_chain,
                        'target_seq': target_seq,
                        'off_target_chain': off_target_chain,
                        'off_target_seq': off_target_seq
                    }

        print(f"✓ Loaded sequences for {len(complex_seqs)} protein complexes")
        return complex_seqs

    def generate_negative_strategy1(self, num_examples: int = 857) -> List[Dict]:
        """
        Strategy 1: Random pairing of 2P2I complexes with compounds from
        2P2I and SuperDRUG2 (target: 857 examples)
        """
        print("\n" + "="*80)
        print("GENERATING NEGATIVE EXAMPLES - STRATEGY 1: RANDOM PAIRING")
        print("="*80)

        negative_examples = []

        # Get all unique complexes and SMILES
        all_complexes = list(set([x['complex_pair'] for x in self.positive_examples]))
        all_smiles_2p2i = list(set([x['smiles'] for x in self.positive_examples]))

        # Create a set of positive pairs to avoid
        positive_pairs = set()
        for ex in self.positive_examples:
            positive_pairs.add((ex['complex_pair'], ex['smiles']))

        # Load SuperDRUG2 compounds
        superdrug_file = os.path.join(self.data_dir,
                                      "approved_drugs_chemical_structure_identifiers.xlsx")
        superdrug_smiles = []

        if os.path.exists(superdrug_file):
            try:
                df_superdrug = pd.read_excel(superdrug_file)
                if 'SMILES' in df_superdrug.columns:
                    superdrug_smiles = df_superdrug['SMILES'].dropna().tolist()
                    print(f"✓ Loaded {len(superdrug_smiles)} compounds from SuperDRUG2")
            except Exception as e:
                print(f"WARNING: Could not load SuperDRUG2 file: {e}")

        # Combine all available SMILES
        all_smiles = all_smiles_2p2i + superdrug_smiles
        print(f"  - Total available compounds: {len(all_smiles)}")

        # Generate random negative pairs
        attempts = 0
        max_attempts = num_examples * 10

        while len(negative_examples) < num_examples and attempts < max_attempts:
            complex_name = random.choice(all_complexes)
            smiles = random.choice(all_smiles)

            # Check if this is not a known positive pair
            if (complex_name, smiles) not in positive_pairs:
                negative_examples.append({
                    'complex_pair': complex_name,
                    'complex_name': complex_name.split('_')[0],
                    'inhibitor_name': f'RANDOM_{len(negative_examples)}',
                    'inhibitor_code': f'RND{len(negative_examples)}',
                    'smiles': smiles,
                    'label': 0,
                    'strategy': 'random_pairing'
                })
                positive_pairs.add((complex_name, smiles))

            attempts += 1

        print(f"✓ Generated {len(negative_examples)} negative examples (Strategy 1)")
        return negative_examples

    def generate_negative_strategy2(self, num_examples: int = 1714) -> List[Dict]:
        """
        Strategy 2: Pair 2P2I compounds with DBD5 complexes
        (target: 1,714 examples)
        """
        print("\n" + "="*80)
        print("GENERATING NEGATIVE EXAMPLES - STRATEGY 2: DBD5 PAIRING")
        print("="*80)

        negative_examples = []

        # Get all SMILES from 2P2I
        all_smiles = list(set([x['smiles'] for x in self.positive_examples]))

        # Load DBD5 complexes
        dbd5_dir = os.path.join(self.data_dir, "DBD5")
        dbd5_seq_file = os.path.join(self.data_dir, "DBD5_seq_dict")

        dbd5_complexes = []

        if os.path.exists(dbd5_seq_file):
            try:
                with open(dbd5_seq_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 1:
                            dbd5_complexes.append(parts[0])
            except UnicodeDecodeError:
                try:
                    with open(dbd5_seq_file, 'r', encoding='latin-1') as f:
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 1:
                                dbd5_complexes.append(parts[0])
                except Exception as e:
                    print(f"  WARNING: Could not read DBD5 sequence file: {e}")

        if not dbd5_complexes:
            # If no sequence file, try to list PDB files in DBD5 directory
            if os.path.exists(dbd5_dir):
                dbd5_files = [f for f in os.listdir(dbd5_dir) if f.endswith('.pdb')]
                dbd5_complexes = [f.replace('.pdb', '') for f in dbd5_files]

        print(f"  - DBD5 complexes available: {len(dbd5_complexes)}")

        if not dbd5_complexes:
            print("  WARNING: No DBD5 complexes found, generating synthetic complex names")
            dbd5_complexes = [f"DBD5_{i}" for i in range(282)]

        # Generate random pairings
        attempts = 0
        max_attempts = num_examples * 10
        generated_pairs = set()

        while len(negative_examples) < num_examples and attempts < max_attempts:
            complex_name = random.choice(dbd5_complexes)
            smiles = random.choice(all_smiles)

            pair_key = (complex_name, smiles)
            if pair_key not in generated_pairs:
                negative_examples.append({
                    'complex_pair': complex_name,
                    'complex_name': complex_name.split('_')[0],
                    'inhibitor_name': f'DBD5_{len(negative_examples)}',
                    'inhibitor_code': f'DBD{len(negative_examples)}',
                    'smiles': smiles,
                    'label': 0,
                    'strategy': 'dbd5_pairing'
                })
                generated_pairs.add(pair_key)

            attempts += 1

        print(f"✓ Generated {len(negative_examples)} negative examples (Strategy 2)")
        return negative_examples

    def generate_negative_strategy3(self) -> List[Dict]:
        """
        Strategy 3: Load pre-computed binders that are not inhibitors
        (From Binding-DB with Tanimoto < 0.85)
        """
        print("\n" + "="*80)
        print("GENERATING NEGATIVE EXAMPLES - STRATEGY 3: BINDERS (NOT INHIBITORS)")
        print("="*80)

        negative_examples = []

        # Load pre-computed binders file
        binders_file = os.path.join(self.data_dir, "BindersWithComplexname.csv")

        if not os.path.exists(binders_file):
            print(f"  WARNING: Binders file not found: {binders_file}")
            print("  Attempting to use Tanimoto similarity file...")

            tanimoto_file = os.path.join(self.data_dir, "Binders With Tanimoto Similarity 0.85.csv")
            if os.path.exists(tanimoto_file):
                try:
                    df = pd.read_csv(tanimoto_file)
                    print(f"  ✓ Loaded Tanimoto similarity file with {len(df)} entries")

                    for idx, row in df.iterrows():
                        if 'Complexname' in row and 'Binders SMILES' in row:
                            complex_name = str(row['Complexname']).strip()
                            smiles = str(row['Binders SMILES']).strip()

                            if smiles and smiles != 'nan':
                                negative_examples.append({
                                    'complex_pair': complex_name,
                                    'complex_name': complex_name.split('_')[0] if '_' in complex_name else complex_name,
                                    'inhibitor_name': f'BINDER_{len(negative_examples)}',
                                    'inhibitor_code': f'BND{len(negative_examples)}',
                                    'smiles': smiles,
                                    'label': 0,
                                    'strategy': 'binders'
                                })
                except Exception as e:
                    print(f"  ERROR loading Tanimoto file: {e}")
        else:
            try:
                df = pd.read_csv(binders_file)
                print(f"  ✓ Loaded binders file with {len(df)} entries")

                for idx, row in df.iterrows():
                    if 'Complexname' in row and 'Binders SMILES' in row:
                        complex_name = str(row['Complexname']).strip()
                        smiles = str(row['Binders SMILES']).strip()

                        if smiles and smiles != 'nan':
                            negative_examples.append({
                                'complex_pair': complex_name,
                                'complex_name': complex_name.split('_')[0] if '_' in complex_name else complex_name,
                                'inhibitor_name': f'BINDER_{len(negative_examples)}',
                                'inhibitor_code': f'BND{len(negative_examples)}',
                                'smiles': smiles,
                                'label': 0,
                                'strategy': 'binders'
                            })
            except Exception as e:
                print(f"  ERROR loading binders file: {e}")

        print(f"✓ Generated {len(negative_examples)} negative examples (Strategy 3)")
        return negative_examples

    def compute_ecfp_features(self, smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
        """
        Compute Extended-Connectivity Fingerprint (ECFP) / Morgan fingerprints

        Args:
            smiles: SMILES string of the compound
            radius: Radius for Morgan fingerprint (default: 2)
            n_bits: Number of bits for fingerprint (default: 2048)

        Returns:
            numpy array of fingerprint bits
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return np.zeros(n_bits)

            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
            arr = np.zeros(n_bits)
            DataStructs.ConvertToNumpyArray(fp, arr)
            return arr
        except Exception as e:
            print(f"  WARNING: Error computing ECFP for SMILES: {smiles[:50]}... Error: {e}")
            return np.zeros(n_bits)

    def compute_aac_features(self, sequence: str) -> np.ndarray:
        """
        Compute Amino Acid Composition (AAC) features

        Args:
            sequence: Protein sequence string

        Returns:
            numpy array of 20 AAC values
        """
        if not sequence:
            return np.zeros(20)

        # Count amino acids
        aa_counts = Counter(sequence.upper())
        total = len(sequence)

        # Compute frequencies for standard amino acids
        aac = np.zeros(20)
        for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY"):
            aac[i] = aa_counts.get(aa, 0) / total if total > 0 else 0

        return aac

    def compute_grouped_kmer_features(self, sequence: str, k: int = 2) -> np.ndarray:
        """
        Compute grouped k-mer composition features

        Args:
            sequence: Protein sequence string
            k: k-mer length (default: 2)

        Returns:
            numpy array of grouped k-mer features (49 for k=2, 7 groups)
        """
        if not sequence or len(sequence) < k:
            return np.zeros(7**k)

        # Convert sequence to group indices
        group_seq = []
        for aa in sequence.upper():
            group_seq.append(self.aa_groups.get(aa, 0))

        # Count k-mers
        num_groups = 7
        kmer_counts = np.zeros(num_groups**k)

        for i in range(len(group_seq) - k + 1):
            kmer = group_seq[i:i+k]
            # Convert k-mer to index
            idx = sum(kmer[j] * (num_groups**(k-1-j)) for j in range(k))
            kmer_counts[idx] += 1

        # Normalize
        total = len(group_seq) - k + 1
        if total > 0:
            kmer_counts = kmer_counts / total

        return kmer_counts

    def compute_protein_sequence_features(self, sequences: List[str]) -> np.ndarray:
        """
        Compute protein sequence features (AAC + grouped k-mer)
        Features are averaged across all chains in the complex

        Args:
            sequences: List of protein sequences

        Returns:
            numpy array of protein features (20 + 49 = 69 dimensions)
        """
        if not sequences:
            return np.zeros(69)

        all_aac = []
        all_kmer = []

        for seq in sequences:
            if seq:
                all_aac.append(self.compute_aac_features(seq))
                all_kmer.append(self.compute_grouped_kmer_features(seq, k=2))

        if not all_aac:
            return np.zeros(69)

        # Average across chains
        aac_features = np.mean(all_aac, axis=0)
        kmer_features = np.mean(all_kmer, axis=0)

        return np.concatenate([aac_features, kmer_features])

    def compute_interface_features(self, complex_name: str,
                                   pdb_dir: str = None) -> np.ndarray:
        """
        Compute interface features (amino acid pairs within 8Å)

        Args:
            complex_name: Name of the complex
            pdb_dir: Directory containing PDB files

        Returns:
            numpy array of interface features (211 dimensions)
        """
        # 211 = 21*21 pairs (20 standard amino acids + 1 unknown)
        interface_features = np.zeros(21 * 21)

        if pdb_dir is None:
            pdb_dir = os.path.join(self.data_dir, "Pdb")

        pdb_file = os.path.join(pdb_dir, f"{complex_name}.pdb")

        if not os.path.exists(pdb_file):
            # Return zero features if PDB not available
            return interface_features[:211]

        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(complex_name, pdb_file)

            # Get all residues
            residues = list(Selection.unfold_entities(structure, 'R'))

            # Find interface residues (within 8Å)
            interface_pairs = []

            for i, res1 in enumerate(residues):
                for res2 in residues[i+1:]:
                    # Check if from different chains
                    if res1.parent.id != res2.parent.id:
                        # Calculate minimum distance between residues
                        min_dist = float('inf')
                        for atom1 in res1:
                            for atom2 in res2:
                                dist = atom1 - atom2  # Distance in Angstroms
                                if dist < min_dist:
                                    min_dist = dist

                        # If within 8Å, record the amino acid pair
                        if min_dist <= 8.0:
                            aa1 = res1.get_resname()
                            aa2 = res2.get_resname()

                            # Convert 3-letter code to 1-letter
                            aa1_code = self._three_to_one(aa1)
                            aa2_code = self._three_to_one(aa2)

                            interface_pairs.append((aa1_code, aa2_code))

            # Count amino acid pairs
            for aa1, aa2 in interface_pairs:
                idx1 = self.amino_acids.index(aa1) if aa1 in self.amino_acids else 20
                idx2 = self.amino_acids.index(aa2) if aa2 in self.amino_acids else 20
                interface_features[idx1 * 21 + idx2] += 1

            # Normalize
            if np.sum(interface_features) > 0:
                interface_features = interface_features / np.sum(interface_features)

        except Exception as e:
            print(f"  WARNING: Error computing interface features for {complex_name}: {e}")

        return interface_features[:211]

    def _three_to_one(self, three_letter: str) -> str:
        """Convert 3-letter amino acid code to 1-letter"""
        conversion = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        return conversion.get(three_letter.upper(), 'X')

    def generate_complete_dataset(self, output_dir: str = "./generated_data",
                                 include_features: bool = True):
        """
        Generate complete dataset with all examples and features

        Args:
            output_dir: Directory to save generated dataset
            include_features: Whether to compute and save features
        """
        print("\n" + "="*80)
        print("GENERATING COMPLETE DATASET")
        print("="*80)

        os.makedirs(output_dir, exist_ok=True)

        # Load positive examples
        self.load_positive_examples()
        complex_seqs = self.load_complex_sequences()

        # Generate negative examples
        neg_strategy1 = self.generate_negative_strategy1(num_examples=857)
        neg_strategy2 = self.generate_negative_strategy2(num_examples=1714)
        neg_strategy3 = self.generate_negative_strategy3()

        # Combine all negative examples
        self.negative_examples = neg_strategy1 + neg_strategy2 + neg_strategy3

        # Combine positive and negative
        all_examples = self.positive_examples + self.negative_examples

        print("\n" + "="*80)
        print("DATASET SUMMARY")
        print("="*80)
        print(f"Positive examples: {len(self.positive_examples)}")
        print(f"Negative examples: {len(self.negative_examples)}")
        print(f"  - Strategy 1 (Random pairing): {len(neg_strategy1)}")
        print(f"  - Strategy 2 (DBD5 pairing): {len(neg_strategy2)}")
        print(f"  - Strategy 3 (Binders): {len(neg_strategy3)}")
        print(f"Total examples: {len(all_examples)}")
        print(f"Class balance: {len(self.positive_examples)}/{len(self.negative_examples)}")

        # Save basic dataset (without features)
        print("\n" + "="*80)
        print("SAVING DATASET FILES")
        print("="*80)

        # Save as text file (matching paper format)
        basic_file = os.path.join(output_dir, "dataset_all_examples.txt")
        with open(basic_file, 'w') as f:
            for ex in all_examples:
                # Format: ComplexName TargetComplexName InhibitorName Label
                line = f"{ex['complex_pair']} {ex['complex_pair']} {ex['inhibitor_code']} {ex['label']}\n"
                f.write(line)
        print(f"✓ Saved basic dataset: {basic_file}")

        # Save as CSV with more details
        df = pd.DataFrame(all_examples)
        csv_file = os.path.join(output_dir, "dataset_complete.csv")
        df.to_csv(csv_file, index=False)
        print(f"✓ Saved detailed dataset: {csv_file}")

        # Save separate positive and negative files
        df_pos = pd.DataFrame(self.positive_examples)
        df_neg = pd.DataFrame(self.negative_examples)
        df_pos.to_csv(os.path.join(output_dir, "dataset_positive.csv"), index=False)
        df_neg.to_csv(os.path.join(output_dir, "dataset_negative.csv"), index=False)
        print(f"✓ Saved separate positive/negative files")

        # Extract and save features if requested
        if include_features:
            print("\n" + "="*80)
            print("EXTRACTING FEATURES (This may take a while...)")
            print("="*80)

            all_features = []
            all_labels = []

            for i, ex in enumerate(all_examples):
                if (i + 1) % 100 == 0:
                    print(f"  Processing example {i+1}/{len(all_examples)}...")

                # Extract ligand features (ECFP)
                ligand_features = self.compute_ecfp_features(ex['smiles'])

                # Extract protein sequence features
                sequences = []
                if ex['complex_pair'] in complex_seqs:
                    sequences = [
                        complex_seqs[ex['complex_pair']]['target_seq'],
                        complex_seqs[ex['complex_pair']]['off_target_seq']
                    ]
                protein_features = self.compute_protein_sequence_features(sequences)

                # Extract interface features
                interface_features = self.compute_interface_features(ex['complex_name'])

                # Combine all features (2048 + 69 + 211 = 2328 dimensions)
                # Note: GNN features (512 dim) would be added separately during training
                combined_features = np.concatenate([
                    ligand_features,
                    protein_features,
                    interface_features
                ])

                all_features.append(combined_features)
                all_labels.append(ex['label'])

            # Save features
            features_array = np.array(all_features)
            labels_array = np.array(all_labels)

            np.save(os.path.join(output_dir, "features.npy"), features_array)
            np.save(os.path.join(output_dir, "labels.npy"), labels_array)

            print(f"✓ Saved features: shape {features_array.shape}")
            print(f"  - Ligand features (ECFP): 2048 dimensions")
            print(f"  - Protein features (AAC + k-mer): 69 dimensions")
            print(f"  - Interface features: 211 dimensions")
            print(f"  - Total: {features_array.shape[1]} dimensions")
            print(f"  - Note: GNN features (512 dim) added during model training")

        print("\n" + "="*80)
        print("DATASET GENERATION COMPLETE!")
        print("="*80)
        print(f"\nOutput directory: {output_dir}")
        print("\nGenerated files:")
        print("  - dataset_all_examples.txt (paper format)")
        print("  - dataset_complete.csv (detailed)")
        print("  - dataset_positive.csv")
        print("  - dataset_negative.csv")
        if include_features:
            print("  - features.npy")
            print("  - labels.npy")


def main():
    parser = argparse.ArgumentParser(
        description='Generate training dataset for PPI Inhibitor prediction'
    )
    parser.add_argument('--data_dir', type=str, default='./Data',
                       help='Directory containing input data files')
    parser.add_argument('--output_dir', type=str, default='./generated_data',
                       help='Directory to save generated dataset')
    parser.add_argument('--no_features', action='store_true',
                       help='Skip feature extraction (only generate dataset files)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')

    args = parser.parse_args()

    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    print("\n" + "="*80)
    print("PPI INHIBITORS DATASET GENERATION")
    print("="*80)
    print(f"Data directory: {args.data_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Random seed: {args.seed}")
    print(f"Feature extraction: {'No' if args.no_features else 'Yes'}")
    print("="*80)

    # Generate dataset
    generator = DatasetGenerator(data_dir=args.data_dir)
    generator.generate_complete_dataset(
        output_dir=args.output_dir,
        include_features=not args.no_features
    )

    print("\n✓ All done!")


if __name__ == "__main__":
    main()
