"""Jaccard index between the miRNAs prioritised by ARM and by SVD.

The two lists below are taken directly from the published paper:
    ARM: Table 2, miRNAs with the strongest association rules (highest lift)
    SVD: Table 3, top 20 miRNAs ranked by SVD score

J(A, S) = |A ∩ S| / |A ∪ S|

Run from the repository root:
    python scripts/jaccard_index.py
"""

ARM_MIRNAS = {
    "miR-141-3p", "miR-200a-3p", "miR-1-3p", "miR-206", "miR-200bc-3p",
    "miR-429", "miR-367-3p", "miR-32-5p", "miR-92-3p", "miR-25-3p",
    "miR-363-3p",
}

SVD_TOP20 = {
    "miR-29-3p", "miR-92-3p", "miR-9-5p", "miR-101-3p.2", "miR-367-3p",
    "miR-32-5p", "miR-219-5p", "miR-373-3p", "miR-302-3p", "miR-372-3p",
    "miR-520-3p", "miR-204-5p", "miR-23-3p", "miR-21-5p", "miR-590-5p",
    "miR-133b", "miR-365-3p", "miR-200a-3p", "miR-141-3p", "miR-133a-3p.2",
}


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b)


def main() -> None:
    shared = sorted(ARM_MIRNAS & SVD_TOP20)
    union = ARM_MIRNAS | SVD_TOP20
    print(f"ARM miRNAs: {len(ARM_MIRNAS)}  SVD miRNAs: {len(SVD_TOP20)}")
    print(f"Shared: {len(shared)}  Union: {len(union)}")
    print(f"Jaccard index: {jaccard(ARM_MIRNAS, SVD_TOP20):.3f}")
    print("Shared miRNAs: " + ", ".join(shared))


if __name__ == "__main__":
    main()
