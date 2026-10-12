import random
from collections import Counter

BASES = "ATCG"


def generate_sequence(length=20):
    """Generate a random DNA sequence of the given length."""
    return "".join(random.choice(BASES) for _ in range(length))


def base_frequency(seq):
    """Count each single letter (A, T, C, G)."""
    counts = Counter(seq)
    return {b: counts.get(b, 0) for b in BASES}


def kmer_frequency(seq, k):
    """Count overlapping combinations of k letters, sorted by frequency (descending)."""
    kmers = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    return sorted(kmers.items(), key=lambda x: (-x[1], x[0]))


def print_table(title, items, total):
    print(f"\n{title}")
    print("-" * 30)
    for item, count in items:
        print(f"{item:<6}{count:>4}   {count / total * 100:6.2f}%")


def main():
    seq = generate_sequence(20)
    print(f"DNA sequence: {seq}")

    singles = sorted(base_frequency(seq).items(), key=lambda x: -x[1])
    print_table("Single letter frequency", singles, len(seq))

    for k in (2, 3):
        result = kmer_frequency(seq, k)
        total = len(seq) - k + 1
        print_table(f"{k}-letter combination frequency (descending)", result, total)


if __name__ == "__main__":
    main()
