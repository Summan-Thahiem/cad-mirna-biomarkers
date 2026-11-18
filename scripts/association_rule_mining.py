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
"""

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


def main() -> None:
    OUT.mkdir(exist_ok=True)
    matrix = build_binary_matrix(DATA)
    print(f"Binary matrix: {matrix.shape[0]} genes x {matrix.shape[1]} miRNAs")

    itemsets = apriori(matrix, min_support=MIN_SUPPORT, use_colnames=True)
    rules = association_rules(itemsets, metric="lift", min_threshold=MIN_LIFT)
    print(f"Frequent itemsets: {len(itemsets)}  Association rules: {len(rules)}")

    itemsets_out = itemsets.assign(
        itemsets=itemsets["itemsets"].apply(joined),
        size=itemsets["itemsets"].apply(len),
    ).sort_values(["size", "support"], ascending=[False, False])
    itemsets_out.to_csv(OUT / "frequent_itemsets.csv", index=False)

    rules_out = rules.assign(
        antecedents=rules["antecedents"].apply(joined),
        consequents=rules["consequents"].apply(joined),
    )[["antecedents", "consequents", "support", "confidence", "lift"]]
    rules_out = rules_out.sort_values("lift", ascending=False)
    rules_out.to_csv(OUT / "association_rules.csv", index=False)

    # One-to-one rules (single miRNA -> single miRNA), strongest first
    pairwise = rules_out[
        ~rules_out["antecedents"].str.contains(",")
        & ~rules_out["consequents"].str.contains(",")
    ]
    pairwise.head(TOP_N_PAIRS).to_csv(OUT / "top_pairwise_rules.csv", index=False)
    print("\nTop pairwise rules by lift:")
    print(pairwise.head(10).to_string(index=False))

    # Heatmap of lift values for the strongest pairwise rules
    top = pairwise.head(TOP_N_PAIRS)
    heat = top.pivot_table(index="antecedents", columns="consequents", values="lift")
    plt.figure(figsize=(11, 9))
    sns.heatmap(heat, annot=True, fmt=".1f", cmap="YlGnBu", linewidths=0.5,
                cbar_kws={"label": "Lift"})
    plt.title("Lift of the strongest pairwise miRNA association rules")
    plt.xlabel("Consequent miRNA")
    plt.ylabel("Antecedent miRNA")
    plt.tight_layout()
    plt.savefig(OUT / "lift_heatmap.png", dpi=200)
    print(f"\nResults written to {OUT}")


if __name__ == "__main__":
    main()
