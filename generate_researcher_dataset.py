"""
Script to generate the researcher's dataset with exact specifications by sampling
from the precomputed dataset WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt:

- TOTAL EXAMPLES: 11,127
- POSITIVE EXAMPLES: 714
- NEGATIVE EXAMPLES (Non-inhibitors): 10,413
- UNIQUE PROTEIN COMPLEXES: 22
- UNIQUE INHIBITORS (from positive examples): 606

This script reads from the precomputed dataset and intelligently samples to match
the exact research specifications while maintaining all unique inhibitors.

Author: Claude Code
Date: 2025-11-12
"""

import random
import os
from collections import defaultdict
from typing import List, Tuple, Set
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
PRECOMPUTED_FILE = os.path.join(DATA_DIR, 'WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt')
OUTPUT_FILE = os.path.join(DATA_DIR, 'researcher_dataset.txt')


def read_precomputed_dataset(filename: str) -> Tuple[List[Tuple], List[Tuple], Set[str], Set[str]]:
    """
    Read the precomputed dataset and separate into positives and negatives.

    Returns: (positive_examples, negative_examples, complexes, positive_inhibitors)
    """
    positive_examples = []
    negative_examples = []
    complexes = set()
    positive_inhibitors = set()

    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    complex_name = parts[0]
                    target_complex = parts[1]
                    # Compound ID can have spaces, so take everything between target_complex and label
                    compound_id = ' '.join(parts[2:-1])
                    label = parts[-1]

                    example = (complex_name, target_complex, compound_id, label)
                    complexes.add(complex_name)

                    if label == '1.0':
                        positive_examples.append(example)
                        positive_inhibitors.add(compound_id)
                    else:
                        negative_examples.append(example)
    except FileNotFoundError:
        print(f"Error: File {filename} not found!")
        sys.exit(1)

    return positive_examples, negative_examples, complexes, positive_inhibitors


def sample_positive_examples(positive_examples: List[Tuple], target_count: int,
                             target_unique_inhibitors: int) -> List[Tuple]:
    """
    Sample positive examples ensuring we maintain all unique inhibitors.

    The precomputed dataset has 857 positives with 606 unique inhibitors.
    We need 714 examples while maintaining all 606 unique inhibitors.
    """
    # Group by inhibitor to ensure we keep at least one example per inhibitor
    inhibitor_examples = defaultdict(list)
    for example in positive_examples:
        inhibitor_id = example[2]
        inhibitor_examples[inhibitor_id].append(example)

    print(f"  Found {len(inhibitor_examples)} unique inhibitors")

    # First pass: select one example for each unique inhibitor
    selected = []
    for inhibitor_id, examples in inhibitor_examples.items():
        # Randomly select one example for this inhibitor
        selected.append(random.choice(examples))

    print(f"  After first pass: {len(selected)} examples (one per inhibitor)")

    # Second pass: add more examples to reach target count
    # Create a pool of remaining examples
    remaining = [ex for ex in positive_examples if ex not in selected]

    if remaining and len(selected) < target_count:
        needed = target_count - len(selected)
        additional = random.sample(remaining, min(needed, len(remaining)))
        selected.extend(additional)

    # Verify we have all unique inhibitors
    selected_inhibitors = set([ex[2] for ex in selected])

    if len(selected_inhibitors) != target_unique_inhibitors:
        print(f"  Warning: Got {len(selected_inhibitors)} unique inhibitors, expected {target_unique_inhibitors}")

    return selected


def sample_negative_examples(negative_examples: List[Tuple], target_count: int) -> List[Tuple]:
    """
    Sample negative examples to reach target count.

    The precomputed dataset has 14,838 negatives. We need 10,413.
    """
    if len(negative_examples) < target_count:
        print(f"  Warning: Only {len(negative_examples)} negatives available, need {target_count}")
        return negative_examples

    return random.sample(negative_examples, target_count)


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
            # Compound ID can have spaces, so take everything between target_complex and label
            compound = ' '.join(parts[2:-1])
            label = parts[-1]

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
    Main function to generate the researcher's dataset from precomputed data.
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

    # Read precomputed dataset
    print("\nStep 1: Reading precomputed dataset...")
    positive_examples, negative_examples, complexes, positive_inhibitors = read_precomputed_dataset(PRECOMPUTED_FILE)
    print(f"  Loaded {len(positive_examples)} positive examples")
    print(f"  Loaded {len(negative_examples)} negative examples")
    print(f"  Found {len(complexes)} unique complexes")
    print(f"  Found {len(positive_inhibitors)} unique inhibitors")

    # Verify source data
    if len(complexes) != TARGET_UNIQUE_COMPLEXES:
        print(f"\n⚠ Warning: Precomputed dataset has {len(complexes)} complexes, expected {TARGET_UNIQUE_COMPLEXES}")
    if len(positive_inhibitors) != TARGET_UNIQUE_INHIBITORS:
        print(f"\n⚠ Warning: Precomputed dataset has {len(positive_inhibitors)} unique inhibitors, expected {TARGET_UNIQUE_INHIBITORS}")

    # Sample positive examples
    print("\nStep 2: Sampling positive examples...")
    sampled_positives = sample_positive_examples(positive_examples, TARGET_POSITIVE,
                                                 TARGET_UNIQUE_INHIBITORS)
    print(f"  Sampled {len(sampled_positives)} positive examples")
    print(f"  Unique inhibitors: {len(set([ex[2] for ex in sampled_positives]))}")
    print(f"  Unique complexes: {len(set([ex[0] for ex in sampled_positives]))}")

    # Sample negative examples
    print("\nStep 3: Sampling negative examples...")
    sampled_negatives = sample_negative_examples(negative_examples, TARGET_NEGATIVE)
    print(f"  Sampled {len(sampled_negatives)} negative examples")

    # Write dataset
    print("\nStep 4: Writing dataset to file...")
    write_dataset(sampled_positives, sampled_negatives, OUTPUT_FILE)

    # Verify dataset
    print("\nStep 5: Verifying dataset...")
    verify_dataset(OUTPUT_FILE)

    print("\n" + "="*60)
    print("DATASET GENERATION COMPLETE")
    print("="*60)
    print(f"\nSource: {PRECOMPUTED_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print("\nYou can now use this dataset for training and evaluation.")


if __name__ == "__main__":
    main()
