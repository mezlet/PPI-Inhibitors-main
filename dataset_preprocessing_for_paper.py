"""
PPI Inhibitors Dataset Preprocessing Script
===========================================

This script preprocesses the PPI inhibitors dataset to match the exact specifications
described in the research paper:

"Predicting small-molecule inhibition of protein complexes"
Yaseen et al., 2024

Dataset Specifications from Paper:
- Training Set: 714 positive examples (inhibitors) from 22 protein complexes
- Negative Examples: ~10,413-14,838 examples from 3 strategies:
  1. Random pairing (2P2I + SuperDRUG2): ~857 examples
  2. 2P2I compounds with DBD5 complexes: ~1,714 examples
  3. Binders that are not inhibitors (BindingDB): ~7,842-11,789 examples
- External Test Set 1: 28 inhibitors from recent publications
- External Test Set 2: 25 SARS-CoV-2 Spike/ACE2 inhibitors

Reference: Table 2 from the paper lists 22 complexes used in LOCO validation
"""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import os

class PPIDatasetPreprocessor:
    """
    Preprocessor for PPI Inhibitors dataset to match research paper specifications.
    """

    # 22 complexes from Table 2 of the research paper (ordered by appearance)
    PAPER_COMPLEXES = [
        '3DAB',    # 1. MDM4/P53
        '3WN7',    # 2. MKEAP1/MNRF2
        '2FLU',    # 3. KEAP1/NRF2
        '1BKD',    # 4. HRAS/SOS1
        '1YCQ',    # 5. XDM2/P53
        '4ESG',    # 6. WDR5/MLL1
        '4QC3',    # 7. BAZ2B/H4
        '3TDU',    # 8. DCN1/UBC12
        '1F47',    # 9. ZIPA/FTSZ
        '2E3K',    # 10. BRD2-2/H4
        '4AJY',    # 11. VHL/HIF1A
        '3D9T',    # 12. CIAP1-BIR3/CASPASE-9
        '2RNY',    # 13. CREBBP/H4
        '3UVW',    # 14. BRD4-1/H4
        '4YY6',    # 15. BRD9/H4
        '1YCR',    # 16. MDM2/P53
        '1BXL',    # 17. BCLXL/BAK
        '2B4J',    # 18. INTEGRASE/LEDGF
        '2XA0',    # 19. BCL2/BAX
        '1Z92',    # 20. IL-2/IL-2R
        '1NW9',    # 21. XIAP-BIR3/SMAC
        '4GQ6',    # 22. MENIN/MLL
    ]

    # External test complexes (not in training)
    EXTERNAL_COMPLEXES = {
        '2dyh': 'MDM2-p53 external',  # External test set 1 (28 inhibitors)
        '6m0j': 'SARS-CoV-2/ACE2'     # External test set 2 (25 inhibitors)
    }

    def __init__(self, data_dir='/home/user/PPI-Inhibitors-main/Data'):
        """
        Initialize the preprocessor.

        Args:
            data_dir: Path to the data directory
        """
        self.data_dir = data_dir
        self.main_dataset_file = os.path.join(data_dir, 'WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt')
        self.external1_file = os.path.join(data_dir, 'External data', '2dyh_all_External_All_Examples.txt')
        self.external2_file = os.path.join(data_dir, 'External data', 'HansonACE2hits_External_All_Examples.txt')

        self.data = None
        self.stats = {}

    def load_main_dataset(self):
        """
        Load the main dataset file.

        Format: Complex_Name Target_Complex Inhibitor_SMILES Label
        """
        print("Loading main dataset...")

        data = []
        with open(self.main_dataset_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split()
                if len(parts) < 3:
                    print(f"Warning: Line {line_num} has fewer than 3 parts: {line}")
                    continue

                # Format: ComplexName TargetComplex InhibitorSMILES Label
                # SMILES can contain spaces, so we need to handle this carefully
                complex_name = parts[0]
                target_complex = parts[1]
                label = float(parts[-1])

                # Everything between target_complex and label is the SMILES
                smiles = ' '.join(parts[2:-1])

                # Extract the base complex ID (e.g., '3UVW' from '3UVW_A_2_B')
                base_complex = complex_name.split('_')[0] if '_' in complex_name else complex_name

                data.append({
                    'complex_name': complex_name,
                    'target_complex': target_complex,
                    'base_complex': base_complex,
                    'smiles': smiles,
                    'label': label,
                    'is_positive': label == 1.0
                })

        self.data = pd.DataFrame(data)
        print(f"Loaded {len(self.data)} examples from main dataset")

        return self.data

    def load_external_datasets(self):
        """
        Load external test datasets.
        """
        print("\nLoading external test datasets...")

        external_data = []

        # Load 2dyh (MDM2-p53) external dataset
        if os.path.exists(self.external1_file):
            with open(self.external1_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    complex_name = parts[0]
                    smiles = ' '.join(parts[1:-1])
                    label = float(parts[-1])

                    external_data.append({
                        'complex_name': complex_name,
                        'target_complex': complex_name,
                        'base_complex': '2dyh',
                        'smiles': smiles,
                        'label': label,
                        'is_positive': label == 1.0,
                        'dataset': 'external_2dyh'
                    })
            print(f"  - Loaded {sum(1 for x in external_data if x['dataset'] == 'external_2dyh')} examples from 2dyh dataset")

        # Load 6m0j (SARS-CoV-2/ACE2) external dataset
        if os.path.exists(self.external2_file):
            start_count = len(external_data)
            with open(self.external2_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    complex_name = parts[0]
                    smiles = ' '.join(parts[1:-1])
                    label = float(parts[-1])

                    external_data.append({
                        'complex_name': complex_name,
                        'target_complex': complex_name,
                        'base_complex': '6m0j',
                        'smiles': smiles,
                        'label': label,
                        'is_positive': label == 1.0,
                        'dataset': 'external_6m0j'
                    })
            print(f"  - Loaded {len(external_data) - start_count} examples from 6m0j dataset")

        return pd.DataFrame(external_data) if external_data else pd.DataFrame()

    def filter_to_paper_complexes(self):
        """
        Filter the dataset to only include the 22 complexes mentioned in the paper.
        """
        print("\nFiltering to paper's 22 complexes...")

        if self.data is None:
            raise ValueError("Data not loaded. Call load_main_dataset() first.")

        # Filter to only the 22 complexes
        filtered_data = self.data[self.data['base_complex'].isin(self.PAPER_COMPLEXES)].copy()

        print(f"  - Original dataset: {len(self.data)} examples")
        print(f"  - Filtered dataset: {len(filtered_data)} examples")
        print(f"  - Removed: {len(self.data) - len(filtered_data)} examples")

        return filtered_data

    def analyze_dataset(self, data=None):
        """
        Analyze the dataset composition and generate statistics.
        """
        if data is None:
            data = self.data

        if data is None:
            raise ValueError("Data not loaded. Call load_main_dataset() first.")

        print("\n" + "="*80)
        print("DATASET ANALYSIS")
        print("="*80)

        # Overall statistics
        total = len(data)
        positives = data['is_positive'].sum()
        negatives = total - positives

        print(f"\nOverall Statistics:")
        print(f"  - Total examples: {total:,}")
        print(f"  - Positive examples (inhibitors): {positives:,} ({100*positives/total:.1f}%)")
        print(f"  - Negative examples: {negatives:,} ({100*negatives/total:.1f}%)")

        # Complex-wise breakdown
        print(f"\nComplex-wise Breakdown:")
        complex_stats = data.groupby('base_complex').agg({
            'is_positive': ['sum', 'count']
        }).round(0)
        complex_stats.columns = ['Positive', 'Total']
        complex_stats['Negative'] = complex_stats['Total'] - complex_stats['Positive']
        complex_stats = complex_stats[['Positive', 'Negative', 'Total']].astype(int)
        complex_stats = complex_stats.sort_values('Total', ascending=False)

        print(complex_stats.to_string())

        # Identify which of the 22 paper complexes are present
        present_complexes = set(data['base_complex'].unique())
        paper_complexes_set = set(self.PAPER_COMPLEXES)

        missing_complexes = paper_complexes_set - present_complexes
        extra_complexes = present_complexes - paper_complexes_set

        if missing_complexes:
            print(f"\n  ⚠ Missing complexes from paper: {missing_complexes}")
        if extra_complexes:
            print(f"\n  ⚠ Extra complexes not in paper: {extra_complexes}")

        # Label distribution
        print(f"\nLabel Distribution:")
        label_counts = data['label'].value_counts().sort_index()
        for label, count in label_counts.items():
            print(f"  - Label {label}: {count:,} examples")

        # Statistics storage
        self.stats = {
            'total': total,
            'positives': positives,
            'negatives': negatives,
            'complex_stats': complex_stats,
            'present_complexes': len(present_complexes),
            'missing_complexes': len(missing_complexes),
            'extra_complexes': len(extra_complexes)
        }

        return self.stats

    def compare_with_paper(self, filtered_data):
        """
        Compare the filtered dataset with the paper's reported statistics.
        """
        print("\n" + "="*80)
        print("COMPARISON WITH RESEARCH PAPER")
        print("="*80)

        paper_stats = {
            'complexes': 22,
            'positives': 714,
            'negatives_total': 14838,  # Based on the actual file
            'negatives_strategy1': 857,
            'negatives_strategy2': 1714,
            'negatives_strategy3': 11789,  # Binders
        }

        actual_positives = filtered_data['is_positive'].sum()
        actual_negatives = len(filtered_data) - actual_positives
        actual_complexes = filtered_data['base_complex'].nunique()

        print(f"\n{'Metric':<40} {'Paper':<15} {'Actual':<15} {'Match':<10}")
        print("-" * 80)
        print(f"{'Number of complexes':<40} {paper_stats['complexes']:<15} {actual_complexes:<15} {'✓' if actual_complexes == paper_stats['complexes'] else '✗'}")
        print(f"{'Positive examples (inhibitors)':<40} {paper_stats['positives']:<15} {actual_positives:<15} {'✓' if actual_positives == paper_stats['positives'] else '✗'}")
        print(f"{'Negative examples':<40} {paper_stats['negatives_total']:<15} {actual_negatives:<15} {'✓' if actual_negatives == paper_stats['negatives_total'] else '✗'}")

        print(f"\nNegative Example Strategy Breakdown (from paper):")
        print(f"  - Strategy 1 (Random pairing): {paper_stats['negatives_strategy1']:,}")
        print(f"  - Strategy 2 (DBD5 complexes): {paper_stats['negatives_strategy2']:,}")
        print(f"  - Strategy 3 (Binders): {paper_stats['negatives_strategy3']:,}")
        print(f"  - Total: {sum([paper_stats['negatives_strategy1'], paper_stats['negatives_strategy2'], paper_stats['negatives_strategy3']]):,}")

        discrepancy = actual_positives - paper_stats['positives']
        if discrepancy != 0:
            print(f"\n⚠ NOTE: Positive examples differ by {abs(discrepancy)} ({'+' if discrepancy > 0 else ''}{discrepancy})")
            if discrepancy > 0:
                print("   The dataset file may contain additional positive examples not used in the paper.")
            else:
                print("   Some positive examples from the paper may be missing from this file.")

        return {
            'match_complexes': actual_complexes == paper_stats['complexes'],
            'match_positives': actual_positives == paper_stats['positives'],
            'match_negatives': actual_negatives == paper_stats['negatives_total']
        }

    def save_preprocessed_data(self, filtered_data, external_data, output_dir=None):
        """
        Save the preprocessed datasets.
        """
        if output_dir is None:
            output_dir = os.path.join(self.data_dir, 'preprocessed')

        os.makedirs(output_dir, exist_ok=True)

        print(f"\nSaving preprocessed data to {output_dir}...")

        # Save training dataset (22 complexes)
        train_file = os.path.join(output_dir, 'training_22_complexes.txt')
        with open(train_file, 'w') as f:
            for _, row in filtered_data.iterrows():
                f.write(f"{row['complex_name']} {row['target_complex']} {row['smiles']} {row['label']}\n")
        print(f"  - Saved training data: {train_file} ({len(filtered_data)} examples)")

        # Save external test sets
        if not external_data.empty:
            # 2dyh dataset
            external1_data = external_data[external_data['dataset'] == 'external_2dyh']
            if not external1_data.empty:
                external1_file = os.path.join(output_dir, 'external_test_2dyh.txt')
                with open(external1_file, 'w') as f:
                    for _, row in external1_data.iterrows():
                        f.write(f"{row['complex_name']} {row['smiles']} {row['label']}\n")
                print(f"  - Saved external test (2dyh): {external1_file} ({len(external1_data)} examples)")

            # 6m0j dataset
            external2_data = external_data[external_data['dataset'] == 'external_6m0j']
            if not external2_data.empty:
                external2_file = os.path.join(output_dir, 'external_test_6m0j.txt')
                with open(external2_file, 'w') as f:
                    for _, row in external2_data.iterrows():
                        f.write(f"{row['complex_name']} {row['smiles']} {row['label']}\n")
                print(f"  - Saved external test (6m0j): {external2_file} ({len(external2_data)} examples)")

        # Save statistics
        stats_file = os.path.join(output_dir, 'dataset_statistics.txt')
        with open(stats_file, 'w') as f:
            f.write("PPI Inhibitors Dataset Statistics\n")
            f.write("="*80 + "\n\n")
            f.write(f"Training Dataset (22 complexes):\n")
            f.write(f"  - Total examples: {len(filtered_data):,}\n")
            f.write(f"  - Positive examples: {filtered_data['is_positive'].sum():,}\n")
            f.write(f"  - Negative examples: {(~filtered_data['is_positive']).sum():,}\n\n")

            if not external_data.empty:
                f.write(f"External Test Datasets:\n")
                for dataset_name in external_data['dataset'].unique():
                    subset = external_data[external_data['dataset'] == dataset_name]
                    f.write(f"  - {dataset_name}: {len(subset):,} examples ")
                    f.write(f"({subset['is_positive'].sum():,} positive, {(~subset['is_positive']).sum():,} negative)\n")

        print(f"  - Saved statistics: {stats_file}")

        return output_dir

    def run_full_preprocessing(self):
        """
        Run the complete preprocessing pipeline.
        """
        print("="*80)
        print("PPI INHIBITORS DATASET PREPROCESSING")
        print("="*80)
        print("\nThis script preprocesses the dataset to match the research paper specifications.")
        print("Reference: 'Predicting small-molecule inhibition of protein complexes'")
        print("           Yaseen et al., 2024")

        # Step 1: Load main dataset
        self.load_main_dataset()

        # Step 2: Analyze original dataset
        print("\n### ORIGINAL DATASET ###")
        self.analyze_dataset()

        # Step 3: Filter to paper's 22 complexes
        filtered_data = self.filter_to_paper_complexes()

        # Step 4: Analyze filtered dataset
        print("\n### FILTERED DATASET (22 complexes) ###")
        self.analyze_dataset(filtered_data)

        # Step 5: Compare with paper
        comparison = self.compare_with_paper(filtered_data)

        # Step 6: Load external datasets
        external_data = self.load_external_datasets()
        if not external_data.empty:
            print("\n### EXTERNAL TEST DATASETS ###")
            self.analyze_dataset(external_data)

        # Step 7: Save preprocessed data
        output_dir = self.save_preprocessed_data(filtered_data, external_data)

        print("\n" + "="*80)
        print("PREPROCESSING COMPLETE")
        print("="*80)
        print(f"\nPreprocessed data saved to: {output_dir}")
        print("\nNext steps:")
        print("  1. Use 'training_22_complexes.txt' for Leave-One-Complex-Out (LOCO) validation")
        print("  2. Use 'external_test_2dyh.txt' for independent validation on MDM2-p53")
        print("  3. Use 'external_test_6m0j.txt' for independent validation on SARS-CoV-2/ACE2")
        print("\n" + "="*80)

        return filtered_data, external_data, output_dir


def main():
    """Main execution function."""
    preprocessor = PPIDatasetPreprocessor()
    filtered_data, external_data, output_dir = preprocessor.run_full_preprocessing()

    return preprocessor, filtered_data, external_data


if __name__ == "__main__":
    preprocessor, filtered_data, external_data = main()
