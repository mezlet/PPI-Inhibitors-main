#!/usr/bin/env python3
"""
Analyze the precomputed PPI Inhibitors dataset and compare to research paper specifications
"""

from collections import Counter

# Load the dataset
dataset_file = '/home/user/PPI-Inhibitors-main/Data/WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt'

print("=" * 80)
print("PPI INHIBITORS DATASET ANALYSIS")
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
                'label': float(label)
            })

print(f"Total examples: {len(data)}")
print()

# Class distribution
print("CLASS DISTRIBUTION:")
print("-" * 80)
pos_count = sum(1 for d in data if d['label'] == 1.0)
neg_count = sum(1 for d in data if d['label'] == 0.0)
print(f"Positive examples (label=1.0): {pos_count}")
print(f"Negative examples (label=0.0): {neg_count}")
print(f"Positive:Negative ratio: 1:{neg_count/pos_count:.2f}")
print()

# Compare to paper
print("COMPARISON TO RESEARCH PAPER:")
print("-" * 80)
print(f"Paper specifications:")
print(f"  - Positive examples: 714")
print(f"  - Negative examples: 10,413")
print(f"  - Total: 11,127")
print(f"  - Pos:Neg ratio: 1:14.6")
print()
print(f"Current dataset:")
print(f"  - Positive examples: {pos_count} (difference: {pos_count - 714:+d})")
print(f"  - Negative examples: {neg_count} (difference: {neg_count - 10413:+d})")
print(f"  - Total: {len(data)} (difference: {len(data) - 11127:+d})")
print(f"  - Pos:Neg ratio: 1:{neg_count/pos_count:.2f}")
print()

# Analyze positive examples
print("POSITIVE EXAMPLES ANALYSIS:")
print("-" * 80)
pos_data = [d for d in data if d['label'] == 1.0]
pos_complexes = [d['complex_name'] for d in pos_data]
unique_complexes_pos = len(set(pos_complexes))
print(f"Unique complexes in positives: {unique_complexes_pos}")
print(f"Paper specification: 22 complexes")
print(f"Difference: {unique_complexes_pos - 22:+d}")
print()

# List complexes with counts
complex_counter = Counter(pos_complexes)
print(f"Top 10 complexes by inhibitor count:")
for complex_name, count in complex_counter.most_common(10):
    print(f"  {complex_name}: {count} inhibitors")
print()

# Analyze negative examples
print("NEGATIVE EXAMPLES ANALYSIS:")
print("-" * 80)
neg_data = [d for d in data if d['label'] == 0.0]

# Try to identify strategies based on complex patterns
# Strategy 1: Random (complex_name == target_complex, inhibitor is numeric or from SuperDRUG2)
# Strategy 2: DBD5 complexes (target_complex different from complex_name, or target_complex not in 2P2I list)
# Strategy 3: Binders (inhibitor_name is numeric, likely from BindingDB)

# Check if complex_name matches target_complex
same_complex_count = sum(1 for d in neg_data if d['complex_name'] == d['target_complex'])
diff_complex_count = sum(1 for d in neg_data if d['complex_name'] != d['target_complex'])

print(f"Negatives where complex_name == target_complex: {same_complex_count}")
print(f"Negatives where complex_name != target_complex: {diff_complex_count}")
print()

# Check if inhibitor names are numeric (likely binders from BindingDB)
numeric_count = sum(1 for d in neg_data if d['inhibitor_name'].isdigit())
non_numeric_count = sum(1 for d in neg_data if not d['inhibitor_name'].isdigit())

print(f"Negatives with numeric inhibitor IDs (likely binders): {numeric_count}")
print(f"Negatives with non-numeric inhibitor IDs: {non_numeric_count}")
print()

# Unique target complexes in negatives
neg_target_complexes = set(d['target_complex'] for d in neg_data)
unique_target_complexes_neg = len(neg_target_complexes)
print(f"Unique target complexes in negatives: {unique_target_complexes_neg}")
print()

# Check for DBD5-style complexes (different pattern)
# Get unique complexes from positives (2P2I complexes)
pos_complexes_set = set(pos_complexes)

# DBD5 complexes would be in negatives but not in positives
dbd5_like_complexes = neg_target_complexes - pos_complexes_set
print(f"Unique complexes in negatives NOT in positives (likely DBD5): {len(dbd5_like_complexes)}")
if len(dbd5_like_complexes) <= 20:
    print(f"Examples: {list(dbd5_like_complexes)[:10]}")
print()

# Estimate strategy breakdown
print("ESTIMATED NEGATIVE STRATEGY BREAKDOWN:")
print("-" * 80)

# Strategy 3: Binders (numeric IDs, same complex)
strategy3_count = sum(1 for d in neg_data if d['complex_name'] == d['target_complex'] and d['inhibitor_name'].isdigit())
print(f"Strategy 3 (Binders - hard negatives): ~{strategy3_count} examples")

# Strategy 2: DBD5 complexes (different target complex)
strategy2_count = sum(1 for d in neg_data if d['complex_name'] != d['target_complex'])
print(f"Strategy 2 (DBD5 complexes): ~{strategy2_count} examples")

# Strategy 1: Random (same complex, non-numeric inhibitor)
strategy1_count = sum(1 for d in neg_data if d['complex_name'] == d['target_complex'] and not d['inhibitor_name'].isdigit())
print(f"Strategy 1 (Random 2P2I/SuperDRUG2): ~{strategy1_count} examples")
print()

# Paper specifications for comparison
print("Paper specifications:")
print(f"  Strategy 1: ~857 examples")
print(f"  Strategy 2: ~1,714 examples")
print(f"  Strategy 3: ~11,789 examples (but only ~8,842 actually used?)")
print()

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("The current precomputed dataset has MORE examples than described in the paper.")
print(f"Current dataset: {len(data)} examples (Pos: {pos_count}, Neg: {neg_count})")
print(f"Paper dataset: 11,127 examples (Pos: 714, Neg: 10,413)")
print()

if pos_count != 714:
    print(f"⚠️  WARNING: Positive examples ({pos_count}) differ from paper (714)")
    print(f"   This might be due to:")
    print(f"   - Different filtering criteria")
    print(f"   - Different version of 2P2I database")
    print(f"   - Including/excluding certain complexes")
    print()

if neg_count > 10413:
    print(f"ℹ️  INFO: More negative examples than paper ({neg_count} vs 10,413)")
    print(f"   This could mean:")
    print(f"   - More negatives were generated (potentially better for training)")
    print(f"   - Different sampling strategy")
    print()

print("To streamline to match paper exactly, we would need to:")
print(f"  1. Filter positives to exactly 714 examples from 22 complexes")
print(f"  2. Adjust negative sampling to get exactly 10,413 negatives")
print(f"  3. Maintain the 1:14.6 ratio")
print()
