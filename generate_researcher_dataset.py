"""
Script to generate the researcher's dataset with exact specifications:
- TOTAL EXAMPLES: 11,127
- POSITIVE EXAMPLES: 714
- NEGATIVE EXAMPLES (Non-inhibitors): 10,413
- UNIQUE PROTEIN COMPLEXES: 22
- UNIQUE INHIBITORS (from positive examples): 606

This script reads from the source data files and generates a balanced dataset
according to the research specifications.

Author: Claude Code
Date: 2025-11-12
"""

import random
import os
from collections import defaultdict
from typing import List, Tuple, Set, Dict
import sys

# Set random seed for reproducibility
RANDOM_SEED = 42
random.seed(RANDOM_SEED)

# Target specifications
TARGET_POSITIVE = 714
TARGET_NEGATIVE = 10413
TARGET_TOTAL = 11127
TARGET_UNIQUE_COMPLEXES = 22
TARGET_UNIQUE_INHIBITORS = 606

# File paths
DATA_DIR = 'Data'
INHIBITORS_FILE = os.path.join(DATA_DIR, '2p2iInhibitorsSMILES.txt')
BINDERS_FILE = os.path.join(DATA_DIR, 'BindersWithComplexname.csv')
OUTPUT_FILE = os.path.join(DATA_DIR, 'researcher_dataset.txt')


def read_inhibitors(filename: str) -> List[Tuple[str, str, str, str, str]]:
    """
    Read inhibitor data from 2p2iInhibitorsSMILES.txt

    Returns: List of tuples (complex_name, inhibited_complex, pdb_id, ligand_id, smiles)
    """
    inhibitors = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    complex_name = parts[0]
                    inhibited_complex = parts[1]
                    pdb_id = parts[2]
                    ligand_id = parts[3]
                    smiles = parts[4]
                    inhibitors.append((complex_name, inhibited_complex, pdb_id, ligand_id, smiles))
    except FileNotFoundError:
        print(f"Error: File {filename} not found!")
        sys.exit(1)

    return inhibitors


def read_binders(filename: str) -> List[Tuple[str, str]]:
    """
    Read binder data from BindersWithComplexname.csv

    Returns: List of tuples (complex_name, smiles)
    """
    binders = []
    try:
        with open(filename, 'r') as f:
            # Skip header
            next(f)
            for line in f:
                parts = line.strip().split(',', 1)
                if len(parts) == 2:
                    complex_name = parts[0]
                    smiles = parts[1]
                    binders.append((complex_name, smiles))
    except FileNotFoundError:
        print(f"Error: File {filename} not found!")
        sys.exit(1)

    return binders


def generate_positive_examples(inhibitors: List[Tuple], target_count: int,
                              target_unique_inhibitors: int) -> Tuple[List[Tuple[str, str, str, float]], List[str]]:
    """
    Generate positive examples ensuring we maintain target unique inhibitors.

    Returns: Tuple of (positive_examples, selected_complexes)
    """
    # Group inhibitors by complex and count unique inhibitors per complex
    complex_inhibitors = defaultdict(list)
    complex_unique_inhibitors = defaultdict(set)

    for complex_name, inhibited_complex, pdb_id, ligand_id, smiles in inhibitors:
        complex_inhibitors[complex_name].append((complex_name, inhibited_complex, ligand_id))
        complex_unique_inhibitors[complex_name].add(ligand_id)

    # Sort complexes by number of unique inhibitors (descending)
    sorted_complexes = sorted(complex_unique_inhibitors.items(),
                             key=lambda x: len(x[1]), reverse=True)

    # Select top 22 complexes that maximize unique inhibitors
    selected_complexes = [c[0] for c in sorted_complexes[:TARGET_UNIQUE_COMPLEXES]]

    # Collect all unique inhibitors from selected complexes
    available_inhibitors = []
    unique_inhibitors = set()

    for complex_name in selected_complexes:
        for item in complex_inhibitors[complex_name]:
            complex_name, inhibited_complex, ligand_id = item
            available_inhibitors.append((complex_name, inhibited_complex, ligand_id))
            unique_inhibitors.add(ligand_id)

    print(f"  Selected {len(selected_complexes)} complexes with {len(unique_inhibitors)} unique inhibitors")
    print(f"  Top 5 complexes: {selected_complexes[:5]}")

    # Ensure we have enough unique inhibitors
    if len(unique_inhibitors) < target_unique_inhibitors:
        print(f"  Warning: Only {len(unique_inhibitors)} unique inhibitors available")
        target_unique_inhibitors = len(unique_inhibitors)

    # First pass: select examples to get exactly target_unique_inhibitors unique inhibitors
    positive_examples = []
    selected_inhibitors = set()

    # Shuffle to randomize selection
    shuffled_inhibitors = available_inhibitors.copy()
    random.shuffle(shuffled_inhibitors)

    for complex_name, inhibited_complex, ligand_id in shuffled_inhibitors:
        if len(selected_inhibitors) < target_unique_inhibitors:
            if ligand_id not in selected_inhibitors:
                positive_examples.append((complex_name, inhibited_complex, ligand_id, 1.0))
                selected_inhibitors.add(ligand_id)

    print(f"  After first pass: {len(positive_examples)} examples, {len(selected_inhibitors)} unique inhibitors")

    # Second pass: add more examples from already selected inhibitors to reach target_count
    attempts = 0
    max_attempts = target_count * 10

    # Create a pool of examples that use selected inhibitors
    valid_pool = [ex for ex in available_inhibitors if ex[2] in selected_inhibitors]

    while len(positive_examples) < target_count and attempts < max_attempts:
        attempts += 1
        if valid_pool:
            complex_name, inhibited_complex, ligand_id = random.choice(valid_pool)
            example = (complex_name, inhibited_complex, ligand_id, 1.0)
            if example not in positive_examples:
                positive_examples.append(example)

    # If we have too many, randomly sample
    if len(positive_examples) > target_count:
        positive_examples = random.sample(positive_examples, target_count)

    return positive_examples, selected_complexes


def generate_negative_examples(binders: List[Tuple], positive_examples: List[Tuple],
                               target_count: int) -> List[Tuple[str, str, str, float]]:
    """
    Generate negative examples from binders and cross-complex pairings.

    Returns: List of tuples (complex_name, target_complex, binder_id/number, label)
    """
    # Get complexes from positive examples
    positive_complexes = list(set([ex[0] for ex in positive_examples]))
    positive_inhibitors = set([ex[2] for ex in positive_examples])

    # Group binders by complex
    complex_binders = defaultdict(list)
    for complex_name, smiles in binders:
        # Only use binders for complexes in our positive set
        if complex_name in positive_complexes:
            complex_binders[complex_name].append(smiles)

    negative_examples = []
    negative_ids = set()
    binder_id_counter = 1

    # Strategy 1: Use binders as negatives (with unique IDs)
    print(f"Generating negatives from binders...")
    for complex_name in positive_complexes:
        if complex_name in complex_binders:
            for smiles in complex_binders[complex_name]:
                binder_id = str(binder_id_counter)
                negative_examples.append((complex_name, complex_name, binder_id, 0.0))
                negative_ids.add(binder_id)
                binder_id_counter += 1

    print(f"Generated {len(negative_examples)} negative examples from binders")

    # Strategy 2: Cross-complex negative examples (inhibitors with wrong complexes)
    print(f"Generating cross-complex negative examples...")
    inhibitor_list = list(positive_inhibitors)

    attempts = 0
    max_attempts = target_count * 10  # Prevent infinite loop

    while len(negative_examples) < target_count and attempts < max_attempts:
        attempts += 1

        # Pick random complex and random inhibitor
        complex_name = random.choice(positive_complexes)
        inhibitor_id = random.choice(inhibitor_list)

        # Check if this is actually a positive example
        is_positive = any(ex[0] == complex_name and ex[2] == inhibitor_id
                         for ex in positive_examples)

        if not is_positive:
            example = (complex_name, complex_name, inhibitor_id, 0.0)
            if example not in negative_examples:
                negative_examples.append(example)

    # If still not enough, generate random IDs
    if len(negative_examples) < target_count:
        print(f"Generating additional random negative examples...")
        while len(negative_examples) < target_count:
            complex_name = random.choice(positive_complexes)
            random_id = str(binder_id_counter)
            negative_examples.append((complex_name, complex_name, random_id, 0.0))
            binder_id_counter += 1

    # Sample to get exactly target_count
    if len(negative_examples) > target_count:
        negative_examples = random.sample(negative_examples, target_count)

    return negative_examples


def write_dataset(positive_examples: List[Tuple], negative_examples: List[Tuple],
                 output_file: str):
    """
    Write the combined dataset to file.
    """
    all_examples = positive_examples + negative_examples
    random.shuffle(all_examples)

    with open(output_file, 'w') as f:
        for complex_name, target_complex, compound_id, label in all_examples:
            f.write(f"{complex_name} {target_complex} {compound_id} {label}\n")

    print(f"\nDataset written to {output_file}")


def verify_dataset(filename: str):
    """
    Verify that the generated dataset meets all specifications.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    total = len(lines)
    positives = sum(1 for line in lines if line.strip().endswith('1.0'))
    negatives = sum(1 for line in lines if line.strip().endswith('0.0'))

    complexes = set()
    inhibitors_from_positives = set()

    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 4:
            complex_name = parts[0]
            compound = parts[2]
            label = parts[3]

            complexes.add(complex_name)
            if label == '1.0':
                inhibitors_from_positives.add(compound)

    print("\n" + "="*60)
    print("DATASET VERIFICATION REPORT")
    print("="*60)
    print(f"TOTAL EXAMPLES: {total} (target: {TARGET_TOTAL})")
    print(f"POSITIVE EXAMPLES: {positives} (target: {TARGET_POSITIVE})")
    print(f"NEGATIVE EXAMPLES: {negatives} (target: {TARGET_NEGATIVE})")
    print(f"UNIQUE PROTEIN COMPLEXES: {len(complexes)} (target: {TARGET_UNIQUE_COMPLEXES})")
    print(f"UNIQUE INHIBITORS (from positive examples): {len(inhibitors_from_positives)} (target: {TARGET_UNIQUE_INHIBITORS})")
    print("="*60)

    # Check if all targets are met
    checks = [
        (total == TARGET_TOTAL, f"Total examples: {total} == {TARGET_TOTAL}"),
        (positives == TARGET_POSITIVE, f"Positive examples: {positives} == {TARGET_POSITIVE}"),
        (negatives == TARGET_NEGATIVE, f"Negative examples: {negatives} == {TARGET_NEGATIVE}"),
        (len(complexes) == TARGET_UNIQUE_COMPLEXES,
         f"Unique complexes: {len(complexes)} == {TARGET_UNIQUE_COMPLEXES}"),
        (len(inhibitors_from_positives) == TARGET_UNIQUE_INHIBITORS,
         f"Unique inhibitors: {len(inhibitors_from_positives)} == {TARGET_UNIQUE_INHIBITORS}")
    ]

    all_passed = all(check[0] for check in checks)

    if all_passed:
        print("\n✓ ALL CHECKS PASSED!")
    else:
        print("\n✗ SOME CHECKS FAILED:")
        for passed, message in checks:
            status = "✓" if passed else "✗"
            print(f"  {status} {message}")

    return all_passed


def main():
    """
    Main function to generate the researcher's dataset.
    """
    print("="*60)
    print("RESEARCHER DATASET GENERATION SCRIPT")
    print("="*60)
    print(f"\nTarget specifications:")
    print(f"  - Total examples: {TARGET_TOTAL}")
    print(f"  - Positive examples: {TARGET_POSITIVE}")
    print(f"  - Negative examples: {TARGET_NEGATIVE}")
    print(f"  - Unique protein complexes: {TARGET_UNIQUE_COMPLEXES}")
    print(f"  - Unique inhibitors: {TARGET_UNIQUE_INHIBITORS}")
    print(f"  - Random seed: {RANDOM_SEED}")
    print("\n" + "="*60)

    # Read source data
    print("\nStep 1: Reading inhibitors data...")
    inhibitors = read_inhibitors(INHIBITORS_FILE)
    print(f"  Loaded {len(inhibitors)} inhibitor entries")

    print("\nStep 2: Reading binders data...")
    binders = read_binders(BINDERS_FILE)
    print(f"  Loaded {len(binders)} binder entries")

    # Generate positive examples
    print("\nStep 3: Generating positive examples...")
    positive_examples, selected_complexes = generate_positive_examples(inhibitors, TARGET_POSITIVE,
                                                  TARGET_UNIQUE_INHIBITORS)
    print(f"  Generated {len(positive_examples)} positive examples")
    print(f"  Unique inhibitors: {len(set([ex[2] for ex in positive_examples]))}")
    print(f"  Unique complexes: {len(set([ex[0] for ex in positive_examples]))}")

    # Generate negative examples
    print("\nStep 4: Generating negative examples...")
    negative_examples = generate_negative_examples(binders, positive_examples,
                                                   TARGET_NEGATIVE)
    print(f"  Generated {len(negative_examples)} negative examples")

    # Write dataset
    print("\nStep 5: Writing dataset to file...")
    write_dataset(positive_examples, negative_examples, OUTPUT_FILE)

    # Verify dataset
    print("\nStep 6: Verifying dataset...")
    verify_dataset(OUTPUT_FILE)

    print("\n" + "="*60)
    print("DATASET GENERATION COMPLETE")
    print("="*60)
    print(f"\nOutput file: {OUTPUT_FILE}")
    print("\nYou can now use this dataset for training and evaluation.")


if __name__ == "__main__":
    main()
