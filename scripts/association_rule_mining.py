"""Association rule mining (Apriori) on the gene-miRNA target matrix.

Builds a binary gene x miRNA matrix from data/mirnas_data_updated_unique.xlsx,
mines frequent miRNA itemsets with the Apriori algorithm (minimum support 0.02)
and generates association rules ranked by lift.

Outputs (written to results/):
    frequent_itemsets.csv
    association_rules.csv
    top_pairwise_rules.csv
    lift_heatmap.png

Run from the repository root:
    python scripts/association_rule_mining.py

Options (defaults match the published analysis):
    --min-support  minimum support for frequent itemsets (default 0.02)
    --min-lift     minimum lift for association rules (default 1.0)
    --top-n        number of pairwise rules to save and plot (default 20)
    --data         path to the gene-miRNA Excel file
    --out          output directory

Example:
    python scripts/association_rule_mining.py --min-support 0.03 --top-n 30
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from mlxtend.frequent_patterns import apriori, association_rules

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "mirnas_data_updated_unique.xlsx"
OUT = ROOT / "results"

MIN_SUPPORT = 0.02
MIN_LIFT = 1.0
TOP_N_PAIRS = 20


def build_binary_matrix(path: Path) -> pd.DataFrame:
    """Return a gene x miRNA matrix with 1 where the miRNA targets the gene."""
    df = pd.read_excel(path)
    df["miRNA"] = df["miRNA"].str.split(",")
    pairs = df.explode("miRNA").reset_index(drop=True)
    # Strip spaces left after commas so "miR-16-5p" and " miR-16-5p" are one miRNA
    pairs["miRNA"] = pairs["miRNA"].str.strip()
    pairs = pairs[pairs["miRNA"] != ""]
    matrix = pd.crosstab(pairs["GENE"], pairs["miRNA"]) > 0
    return matrix


def joined(items) -> str:
    return ", ".join(sorted(items))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Mine frequent miRNA itemsets and association rules from gene-miRNA target data."
    )
    parser.add_argument("--min-support", type=float, default=MIN_SUPPORT,
                        help=f"minimum support for frequent itemsets (default {MIN_SUPPORT})")
    parser.add_argument("--min-lift", type=float, default=MIN_LIFT,
                        help=f"minimum lift for association rules (default {MIN_LIFT})")
    parser.add_argument("--top-n", type=int, default=TOP_N_PAIRS,
                        help=f"number of pairwise rules to save and plot (default {TOP_N_PAIRS})")
    parser.add_argument("--data", type=Path, default=DATA,
                        help="path to the gene-miRNA Excel file")
    parser.add_argument("--out", type=Path, default=OUT,
                        help="output directory")
    args = parser.parse_args()
    if not 0 < args.min_support <= 1:
        parser.error("--min-support must be between 0 and 1")
    if args.top_n < 1:
        parser.error("--top-n must be at least 1")
    return args


def main() -> None:
    args = parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    matrix = build_binary_matrix(args.data)
    print(f"Binary matrix: {matrix.shape[0]} genes x {matrix.shape[1]} miRNAs")
    print(f"Settings: min support {args.min_support}, min lift {args.min_lift}, top {args.top_n} pairs")

    itemsets = apriori(matrix, min_support=args.min_support, use_colnames=True)
    rules = association_rules(itemsets, metric="lift", min_threshold=args.min_lift)
    print(f"Frequent itemsets: {len(itemsets)}  Association rules: {len(rules)}")

    itemsets_out = itemsets.assign(
        itemsets=itemsets["itemsets"].apply(joined),
        size=itemsets["itemsets"].apply(len),
    ).sort_values(["size", "support"], ascending=[False, False])
    itemsets_out.to_csv(out / "frequent_itemsets.csv", index=False)

    rules_out = rules.assign(
        antecedents=rules["antecedents"].apply(joined),
        consequents=rules["consequents"].apply(joined),
    )[["antecedents", "consequents", "support", "confidence", "lift"]]
    rules_out = rules_out.sort_values("lift", ascending=False)
    rules_out.to_csv(out / "association_rules.csv", index=False)

    # One-to-one rules (single miRNA -> single miRNA), strongest first
    pairwise = rules_out[
        ~rules_out["antecedents"].str.contains(",")
        & ~rules_out["consequents"].str.contains(",")
    ]
    pairwise.head(args.top_n).to_csv(out / "top_pairwise_rules.csv", index=False)
    print("\nTop pairwise rules by lift:")
    print(pairwise.head(10).to_string(index=False))

    # Heatmap of lift values for the strongest pairwise rules
    top = pairwise.head(args.top_n)
    heat = top.pivot_table(index="antecedents", columns="consequents", values="lift")
    plt.figure(figsize=(11, 9))
    sns.heatmap(heat, annot=True, fmt=".1f", cmap="YlGnBu", linewidths=0.5,
                cbar_kws={"label": "Lift"})
    plt.title("Lift of the strongest pairwise miRNA association rules")
    plt.xlabel("Consequent miRNA")
    plt.ylabel("Antecedent miRNA")
    plt.tight_layout()
    plt.savefig(out / "lift_heatmap.png", dpi=200)
    print(f"\nResults written to {out}")


if __name__ == "__main__":
    main()
