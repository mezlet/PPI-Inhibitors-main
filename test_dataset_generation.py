"""
Test script for dataset generation
This script performs a quick test of the dataset generation functionality
without computing full features (which can be time-consuming).
"""

import os
import sys
import pandas as pd

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from generate_training_dataset import DatasetGenerator


def test_dataset_generation():
    """Test the dataset generation pipeline"""

    print("="*80)
    print("TESTING DATASET GENERATION")
    print("="*80)

    # Initialize generator
    data_dir = "./Data"
    output_dir = "./test_generated_data"

    if not os.path.exists(data_dir):
        print(f"ERROR: Data directory not found: {data_dir}")
        return False

    generator = DatasetGenerator(data_dir=data_dir)

    # Test 1: Load positive examples
    print("\n[TEST 1] Loading positive examples...")
    try:
        positive_examples = generator.load_positive_examples()
        print(f"✓ Successfully loaded {len(positive_examples)} positive examples")

        if len(positive_examples) == 0:
            print("✗ ERROR: No positive examples loaded")
            return False

        # Check expected range
        if len(positive_examples) < 700 or len(positive_examples) > 750:
            print(f"⚠ WARNING: Expected ~714 positive examples, got {len(positive_examples)}")
        else:
            print(f"✓ Example count in expected range")

    except Exception as e:
        print(f"✗ ERROR loading positive examples: {e}")
        return False

    # Test 2: Load complex sequences
    print("\n[TEST 2] Loading complex sequences...")
    try:
        complex_seqs = generator.load_complex_sequences()
        print(f"✓ Successfully loaded sequences for {len(complex_seqs)} complexes")

        if len(complex_seqs) == 0:
            print("⚠ WARNING: No complex sequences loaded")

    except Exception as e:
        print(f"✗ ERROR loading complex sequences: {e}")
        return False

    # Test 3: Generate negative examples (Strategy 1)
    print("\n[TEST 3] Generating negative examples - Strategy 1 (small sample)...")
    try:
        neg_strategy1 = generator.generate_negative_strategy1(num_examples=50)
        print(f"✓ Successfully generated {len(neg_strategy1)} negative examples")

        if len(neg_strategy1) < 45:
            print("⚠ WARNING: Generated fewer than expected examples")
        else:
            print(f"✓ Generation successful")

    except Exception as e:
        print(f"✗ ERROR generating Strategy 1 negatives: {e}")
        return False

    # Test 4: Generate negative examples (Strategy 2)
    print("\n[TEST 4] Generating negative examples - Strategy 2 (small sample)...")
    try:
        neg_strategy2 = generator.generate_negative_strategy2(num_examples=50)
        print(f"✓ Successfully generated {len(neg_strategy2)} negative examples")

        if len(neg_strategy2) < 45:
            print("⚠ WARNING: Generated fewer than expected examples")
        else:
            print(f"✓ Generation successful")

    except Exception as e:
        print(f"✗ ERROR generating Strategy 2 negatives: {e}")
        return False

    # Test 5: Generate negative examples (Strategy 3)
    print("\n[TEST 5] Generating negative examples - Strategy 3...")
    try:
        neg_strategy3 = generator.generate_negative_strategy3()
        print(f"✓ Successfully generated {len(neg_strategy3)} negative examples")

        if len(neg_strategy3) == 0:
            print("⚠ WARNING: No Strategy 3 negatives generated (binders file may be missing)")
        else:
            print(f"✓ Generation successful")

    except Exception as e:
        print(f"✗ ERROR generating Strategy 3 negatives: {e}")
        return False

    # Test 6: Feature extraction (single example)
    print("\n[TEST 6] Testing feature extraction...")
    try:
        # Test ECFP
        test_smiles = positive_examples[0]['smiles']
        ecfp_features = generator.compute_ecfp_features(test_smiles)
        print(f"✓ ECFP features: {ecfp_features.shape}")

        if ecfp_features.shape[0] != 2048:
            print(f"✗ ERROR: Expected 2048 ECFP features, got {ecfp_features.shape[0]}")
            return False

        # Test AAC
        test_seq = "ACDEFGHIKLMNPQRSTVWY"
        aac_features = generator.compute_aac_features(test_seq)
        print(f"✓ AAC features: {aac_features.shape}")

        if aac_features.shape[0] != 20:
            print(f"✗ ERROR: Expected 20 AAC features, got {aac_features.shape[0]}")
            return False

        # Test k-mer
        kmer_features = generator.compute_grouped_kmer_features(test_seq, k=2)
        print(f"✓ K-mer features: {kmer_features.shape}")

        if kmer_features.shape[0] != 49:
            print(f"✗ ERROR: Expected 49 k-mer features, got {kmer_features.shape[0]}")
            return False

        # Test protein sequence features
        protein_features = generator.compute_protein_sequence_features([test_seq])
        print(f"✓ Protein features: {protein_features.shape}")

        if protein_features.shape[0] != 69:
            print(f"✗ ERROR: Expected 69 protein features, got {protein_features.shape[0]}")
            return False

        print("✓ All feature extraction tests passed")

    except Exception as e:
        print(f"✗ ERROR in feature extraction: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Test 7: Save dataset files
    print("\n[TEST 7] Saving test dataset...")
    try:
        os.makedirs(output_dir, exist_ok=True)

        # Combine small sample
        generator.positive_examples = positive_examples[:50]  # Small sample
        generator.negative_examples = (neg_strategy1[:50] +
                                      neg_strategy2[:50] +
                                      neg_strategy3[:50])

        all_examples = generator.positive_examples + generator.negative_examples

        # Save as CSV
        df = pd.DataFrame(all_examples)
        test_file = os.path.join(output_dir, "test_dataset.csv")
        df.to_csv(test_file, index=False)
        print(f"✓ Saved test dataset: {test_file}")
        print(f"  - Total examples: {len(df)}")
        print(f"  - Positive: {len(df[df['label']==1])}")
        print(f"  - Negative: {len(df[df['label']==0])}")

        # Save in paper format
        text_file = os.path.join(output_dir, "test_dataset.txt")
        with open(text_file, 'w') as f:
            for ex in all_examples:
                line = f"{ex['complex_pair']} {ex['complex_pair']} {ex['inhibitor_code']} {ex['label']}\n"
                f.write(line)
        print(f"✓ Saved text format: {text_file}")

    except Exception as e:
        print(f"✗ ERROR saving dataset: {e}")
        return False

    # All tests passed
    print("\n" + "="*80)
    print("✓ ALL TESTS PASSED!")
    print("="*80)
    print(f"\nTest output saved to: {output_dir}")
    print("\nYou can now run the full dataset generation:")
    print("  python generate_training_dataset.py")
    print("\n")

    return True


if __name__ == "__main__":
    success = test_dataset_generation()
    sys.exit(0 if success else 1)
