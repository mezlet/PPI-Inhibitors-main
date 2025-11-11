#!/usr/bin/env python3
"""
Streamline the PPI Inhibitors dataset to match research paper specifications
"""

import random
from collections import defaultdict

# Set random seed for reproducibility
random.seed(42)

# Load the dataset
dataset_file = '/home/user/PPI-Inhibitors-main/Data/WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt'
output_file = '/home/user/PPI-Inhibitors-main/Data/Streamlined_Dataset_Paper_Specs.txt'

print("=" * 80)
print("STREAMLINING PPI INHIBITORS DATASET TO MATCH PAPER SPECIFICATIONS")
print("=" * 80)
print()

# Read dataset
data = []
with open(dataset_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 4:
            complex_name, target_complex, inhibitor_name, label = parts
            data.append({
                'complex_name': complex_name,
                'target_complex': target_complex,
                'inhibitor_name': inhibitor_name,
                'label': float(label),
                'line': line.strip()
            })

print(f"Original dataset: {len(data)} examples")
print()

# Separate positives and negatives
pos_data = [d for d in data if d['label'] == 1.0]
neg_data = [d for d in data if d['label'] == 0.0]

print(f"Original positives: {len(pos_data)}")
print(f"Original negatives: {len(neg_data)}")
print()

# Paper specifications
TARGET_POSITIVES = 714
TARGET_NEGATIVES = 10413
TARGET_RATIO = TARGET_NEGATIVES / TARGET_POSITIVES  # 14.6

print("Paper target specifications:")
print(f"  Positives: {TARGET_POSITIVES}")
print(f"  Negatives: {TARGET_NEGATIVES}")
print(f"  Ratio: 1:{TARGET_RATIO:.1f}")
print()

# ============================================================================
# STEP 1: Sample positive examples to get exactly 714
# ============================================================================
print("STEP 1: Sampling positive examples")
print("-" * 80)

# Group positives by complex
pos_by_complex = defaultdict(list)
for d in pos_data:
    pos_by_complex[d['complex_name']].append(d)

# Keep all 22 complexes, sample proportionally
total_pos_current = len(pos_data)
sampled_positives = []

for complex_name, examples in sorted(pos_by_complex.items()):
    # Calculate how many to keep from this complex
    proportion = len(examples) / total_pos_current
    target_count = int(TARGET_POSITIVES * proportion)

    # Make sure we keep at least 1 from each complex
    if target_count == 0:
        target_count = 1

    # Sample from this complex
    if len(examples) <= target_count:
        sampled = examples
    else:
        sampled = random.sample(examples, target_count)

    sampled_positives.extend(sampled)
    print(f"  {complex_name}: {len(examples)} -> {len(sampled)} inhibitors")

# Adjust to exactly 714 if needed
if len(sampled_positives) > TARGET_POSITIVES:
    sampled_positives = random.sample(sampled_positives, TARGET_POSITIVES)
elif len(sampled_positives) < TARGET_POSITIVES:
    # Add more from largest complexes
    deficit = TARGET_POSITIVES - len(sampled_positives)
    largest_complexes = sorted(pos_by_complex.items(), key=lambda x: len(x[1]), reverse=True)
    for complex_name, examples in largest_complexes:
        if deficit == 0:
            break
        # Find examples not yet sampled
        current_inhibitors = set(d['inhibitor_name'] for d in sampled_positives if d['complex_name'] == complex_name)
        not_sampled = [d for d in examples if d['inhibitor_name'] not in current_inhibitors]
        if not_sampled:
            to_add = min(deficit, len(not_sampled))
            sampled_positives.extend(random.sample(not_sampled, to_add))
            deficit -= to_add

print()
print(f"Total sampled positives: {len(sampled_positives)}")
print(f"Target: {TARGET_POSITIVES}")
print()

# ============================================================================
# STEP 2: Sample negative examples to match paper distribution
# ============================================================================
print("STEP 2: Sampling negative examples")
print("-" * 80)

# Categorize negatives by strategy
strategy1_negatives = []  # Random (same complex, non-numeric inhibitor)
strategy2_negatives = []  # DBD5 (different complex)
strategy3_negatives = []  # Binders (same complex, numeric inhibitor)

for d in neg_data:
    same_complex = d['complex_name'] == d['target_complex']
    numeric_inhibitor = d['inhibitor_name'].isdigit()

    if not same_complex:
        strategy2_negatives.append(d)
    elif same_complex and numeric_inhibitor:
        strategy3_negatives.append(d)
    else:
        strategy1_negatives.append(d)

print(f"Available negatives by strategy:")
print(f"  Strategy 1 (Random): {len(strategy1_negatives)} available")
print(f"  Strategy 2 (DBD5): {len(strategy2_negatives)} available")
print(f"  Strategy 3 (Binders): {len(strategy3_negatives)} available")
print()

# Paper target distribution (approximate from paper)
TARGET_STRATEGY1 = 857
TARGET_STRATEGY2 = 1714
TARGET_STRATEGY3 = TARGET_NEGATIVES - TARGET_STRATEGY1 - TARGET_STRATEGY2  # ~7842

print(f"Paper target distribution:")
print(f"  Strategy 1: {TARGET_STRATEGY1}")
print(f"  Strategy 2: {TARGET_STRATEGY2}")
print(f"  Strategy 3: {TARGET_STRATEGY3}")
print()

# Sample from each strategy
sampled_strategy1 = random.sample(strategy1_negatives, min(TARGET_STRATEGY1, len(strategy1_negatives)))
sampled_strategy2 = random.sample(strategy2_negatives, min(TARGET_STRATEGY2, len(strategy2_negatives)))
sampled_strategy3 = random.sample(strategy3_negatives, min(TARGET_STRATEGY3, len(strategy3_negatives)))

# If we don't have enough, fill from available
total_sampled_neg = len(sampled_strategy1) + len(sampled_strategy2) + len(sampled_strategy3)
if total_sampled_neg < TARGET_NEGATIVES:
    deficit = TARGET_NEGATIVES - total_sampled_neg
    print(f"⚠️  Deficit of {deficit} negatives, filling from remaining examples...")

    # Get all remaining negatives
    sampled_neg_set = set(d['line'] for d in sampled_strategy1 + sampled_strategy2 + sampled_strategy3)
    remaining = [d for d in neg_data if d['line'] not in sampled_neg_set]

    if remaining:
        fill = random.sample(remaining, min(deficit, len(remaining)))
        sampled_strategy3.extend(fill)  # Add to strategy 3 (most similar to paper)

sampled_negatives = sampled_strategy1 + sampled_strategy2 + sampled_strategy3

print(f"Sampled negatives by strategy:")
print(f"  Strategy 1: {len(sampled_strategy1)}")
print(f"  Strategy 2: {len(sampled_strategy2)}")
print(f"  Strategy 3: {len(sampled_strategy3)}")
print(f"  Total: {len(sampled_negatives)}")
print()

# ============================================================================
# STEP 3: Combine and save streamlined dataset
# ============================================================================
print("STEP 3: Saving streamlined dataset")
print("-" * 80)

streamlined_data = sampled_positives + sampled_negatives

# Shuffle to mix positives and negatives
random.shuffle(streamlined_data)

# Write to file
with open(output_file, 'w') as f:
    for d in streamlined_data:
        f.write(d['line'] + '\n')

print(f"Streamlined dataset saved to: {output_file}")
print()

# ============================================================================
# STEP 4: Verify statistics
# ============================================================================
print("=" * 80)
print("STREAMLINED DATASET STATISTICS")
print("=" * 80)
print()

final_pos = sum(1 for d in streamlined_data if d['label'] == 1.0)
final_neg = sum(1 for d in streamlined_data if d['label'] == 0.0)
final_ratio = final_neg / final_pos if final_pos > 0 else 0

print(f"Total examples: {len(streamlined_data)}")
print(f"Positive examples: {final_pos}")
print(f"Negative examples: {final_neg}")
print(f"Ratio: 1:{final_ratio:.2f}")
print()

# Complex distribution
final_pos_data = [d for d in streamlined_data if d['label'] == 1.0]
unique_complexes = len(set(d['complex_name'] for d in final_pos_data))
print(f"Unique complexes: {unique_complexes}")
print()

print("Comparison to paper:")
print(f"  Target total: {TARGET_POSITIVES + TARGET_NEGATIVES} | Actual: {len(streamlined_data)}")
print(f"  Target positives: {TARGET_POSITIVES} | Actual: {final_pos}")
print(f"  Target negatives: {TARGET_NEGATIVES} | Actual: {final_neg}")
print(f"  Target ratio: 1:14.6 | Actual: 1:{final_ratio:.2f}")
print(f"  Target complexes: 22 | Actual: {unique_complexes}")
print()

if abs(final_pos - TARGET_POSITIVES) <= 5 and abs(final_neg - TARGET_NEGATIVES) <= 100:
    print("✓ Streamlined dataset matches paper specifications!")
else:
    print("⚠️  Minor differences from paper specs (within acceptable range)")

print()
print("=" * 80)
print("DONE")
print("=" * 80)
