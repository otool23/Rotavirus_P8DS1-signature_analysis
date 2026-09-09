#!/usr/bin/env python3

import csv
import math
from pathlib import Path

from scipy.stats import fisher_exact


# ============================================================
# HARD-CODED INPUT / OUTPUT PATHS
# ============================================================

BASE_DIR = Path("/deac/bio/esstmanGrp/otool23/Fisher_test")

VP4_FILE = BASE_DIR / "VP4_fisher_test.csv"
NON_VP4_FILE = BASE_DIR / "VP1_NSP4_fisher_test.csv"

OUTPUT_FILE = BASE_DIR / "Table_S4_Fisher_allele_enrichment.csv"

ALPHA = 0.05

GENE_ORDER = [
    "VP4",
    "VP6",
    "VP1",
    "VP2",
    "VP3",
    "NSP1",
    "NSP2",
    "NSP3",
    "NSP4",
]

VALID_ALLELES = {"A", "B", "C", "D"}


# ============================================================
# HELPERS
# ============================================================

def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError(f"No header row found in {path}")

        rows = []
        for row in reader:
            clean = {
                str(k).strip(): (str(v).strip() if v is not None else "")
                for k, v in row.items()
            }
            rows.append(clean)

    return rows


def normalize_allele(value):
    allele = str(value).strip().upper()
    return allele if allele in VALID_ALLELES else None


def allele_counts(rows, group, gene):
    """
    Count only explicitly assigned A/B/C/D alleles.

    Blank, uncolored, missing, or otherwise unassigned cells are ignored
    and are NOT included in the denominator.
    """
    counts = {"A": 0, "B": 0, "C": 0, "D": 0}

    for row in rows:
        if row.get("Group", "").strip() != group:
            continue

        allele = normalize_allele(row.get(gene, ""))

        if allele is not None:
            counts[allele] += 1

    return counts


def observed_non_b_alleles(test_counts, control_counts):
    """
    Return only non-B allele classes actually observed for that segment.

    Examples:
        VP6  -> A
        VP1  -> A+C+D
        VP3  -> A+C
    """
    observed = []

    for allele in ["A", "C", "D"]:
        if test_counts[allele] > 0 or control_counts[allele] > 0:
            observed.append(allele)

    return observed


def bh_fdr(p_values):
    """
    Benjamini-Hochberg false-discovery-rate correction
    across all 9 gene-level Fisher tests.
    """
    m = len(p_values)

    ranked = sorted(
        enumerate(p_values),
        key=lambda item: item[1]
    )

    adjusted_ranked = [0.0] * m

    for rank, (_, p_value) in enumerate(ranked, start=1):
        adjusted_ranked[rank - 1] = min(
            1.0,
            p_value * m / rank
        )

    # Enforce monotonicity
    for i in range(m - 2, -1, -1):
        adjusted_ranked[i] = min(
            adjusted_ranked[i],
            adjusted_ranked[i + 1]
        )

    adjusted = [0.0] * m

    for (original_index, _), adjusted_p in zip(ranked, adjusted_ranked):
        adjusted[original_index] = adjusted_p

    return adjusted


def n_over_N_percent(n, N):
    """
    Manuscript-ready explicit denominator:
        n/N (percent)
    """
    if N == 0:
        return ""
    return f"{n}/{N} ({100.0 * n / N:.1f}%)"


def format_p(value):
    """
    Keep full numerical precision in the CSV while making tiny values readable.
    """
    if value is None or math.isnan(value):
        return ""
    if value < 0.001:
        return f"{value:.3e}"
    return f"{value:.6f}"


def format_odds_ratio(value):
    if math.isinf(value):
        return "inf"
    if math.isnan(value):
        return ""
    return f"{value:.2f}"


# ============================================================
# VP4 ANALYSIS
# ============================================================

def analyze_vp4(rows):
    """
    VP4 is analyzed separately.

    Biological comparison:
        P[8]-B versus P[8]-A

    Populations:
        P[8]DS-1 versus P[8]Wa

    Only explicitly assigned P[8]-A or P[8]-B alleles are included.
    Uncolored/unassigned VP4 strains are excluded from the denominator.
    """
    test_group = "P[8]DS-1"
    control_group = "P[8]Wa"

    test_counts = allele_counts(rows, test_group, "VP4")
    control_counts = allele_counts(rows, control_group, "VP4")

    test_A = test_counts["A"]
    test_B = test_counts["B"]

    control_A = control_counts["A"]
    control_B = control_counts["B"]

    test_total = test_A + test_B
    control_total = control_A + control_B

    if test_total == 0:
        raise ValueError(
            "VP4: no assigned P[8]-A/P[8]-B alleles found for P[8]DS-1."
        )

    if control_total == 0:
        raise ValueError(
            "VP4: no assigned P[8]-A/P[8]-B alleles found for P[8]Wa."
        )

    odds_ratio, raw_p = fisher_exact(
        [
            [test_B, test_A],
            [control_B, control_A],
        ],
        alternative="two-sided"
    )

    return {
        "Segment": "VP4",
        "Allele_comparison": "P[8]-B vs P[8]-A",
        "Test_group": test_group,
        "P8DS1_n": test_B,
        "P8DS1_N": test_total,
        "Control_group": control_group,
        "Control_n": control_B,
        "Control_N": control_total,
        "Odds_ratio": float(odds_ratio),
        "Fisher_P_raw": float(raw_p),
    }


# ============================================================
# NON-VP4 ANALYSIS
# ============================================================

def analyze_non_vp4_gene(rows, gene):
    """
    For VP6, VP1-VP3, and NSP1-NSP4:

    Test group:
        P[8]DS-1

    Control group:
        P[4]DS-1

    Fisher comparison:
        type-B versus the non-B allele classes actually present
        for that segment.

    Examples:
        A/B only       -> B vs A
        A/B/C          -> B vs A+C
        A/B/C/D        -> B vs A+C+D
    """
    test_group = "P[8]DS-1"
    control_group = "P[4]DS-1"

    test_counts = allele_counts(rows, test_group, gene)
    control_counts = allele_counts(rows, control_group, gene)

    test_total = sum(test_counts.values())
    control_total = sum(control_counts.values())

    if test_total == 0:
        raise ValueError(f"{gene}: no assigned alleles found for {test_group}.")

    if control_total == 0:
        raise ValueError(f"{gene}: no assigned alleles found for {control_group}.")

    non_b_alleles = observed_non_b_alleles(
        test_counts,
        control_counts
    )

    if not non_b_alleles:
        raise ValueError(
            f"{gene}: no non-B comparison allele was observed."
        )

    test_B = test_counts["B"]
    control_B = control_counts["B"]

    test_nonB = sum(
        test_counts[allele]
        for allele in non_b_alleles
    )

    control_nonB = sum(
        control_counts[allele]
        for allele in non_b_alleles
    )

    odds_ratio, raw_p = fisher_exact(
        [
            [test_B, test_nonB],
            [control_B, control_nonB],
        ],
        alternative="two-sided"
    )

    comparison_text = "+".join(non_b_alleles)

    return {
        "Segment": gene,
        "Allele_comparison": f"type-B vs {comparison_text}",
        "Test_group": test_group,
        "P8DS1_n": test_B,
        "P8DS1_N": test_total,
        "Control_group": control_group,
        "Control_n": control_B,
        "Control_N": control_total,
        "Odds_ratio": float(odds_ratio),
        "Fisher_P_raw": float(raw_p),
    }


# ============================================================
# MAIN
# ============================================================

def main():

    vp4_rows = read_csv(VP4_FILE)
    non_vp4_rows = read_csv(NON_VP4_FILE)

    results = []

    # VP4: P[8]-B versus P[8]-A
    results.append(
        analyze_vp4(vp4_rows)
    )

    # All other genes: type-B versus observed non-B allele classes
    for gene in [
        "VP6",
        "VP1",
        "VP2",
        "VP3",
        "NSP1",
        "NSP2",
        "NSP3",
        "NSP4",
    ]:
        results.append(
            analyze_non_vp4_gene(
                non_vp4_rows,
                gene
            )
        )

    # Benjamini-Hochberg correction across all 9 gene-level tests
    raw_p_values = [
        row["Fisher_P_raw"]
        for row in results
    ]

    adjusted_p_values = bh_fdr(
        raw_p_values
    )

    for row, adjusted_p in zip(results, adjusted_p_values):

        row["BH_adjusted_P"] = adjusted_p

        row["Significant_raw"] = (
            "YES"
            if row["Fisher_P_raw"] < ALPHA
            else "NO"
        )

        row["Significant_after_BH"] = (
            "YES"
            if adjusted_p < ALPHA
            else "NO"
        )

        # Manuscript-ready explicit denominator columns
        row["P8DS1_n_over_N_percent"] = n_over_N_percent(
            row["P8DS1_n"],
            row["P8DS1_N"]
        )

        row["Control_n_over_N_percent"] = n_over_N_percent(
            row["Control_n"],
            row["Control_N"]
        )

    order = {
        gene: i
        for i, gene in enumerate(GENE_ORDER)
    }

    results.sort(
        key=lambda row: order[row["Segment"]]
    )

    # ========================================================
    # MANUSCRIPT-READY TABLE S4 OUTPUT
    # ========================================================

    fields = [
        "Segment",
        "Allele_comparison",
        "P8DS1_n_over_N_percent",
        "Control_group",
        "Control_n_over_N_percent",
        "Odds_ratio",
        "Fisher_P_raw",
        "BH_adjusted_P",
        "Significant_raw",
        "Significant_after_BH",
    ]

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with OUTPUT_FILE.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields
        )

        writer.writeheader()

        for row in results:
            writer.writerow({
                "Segment": row["Segment"],
                "Allele_comparison": row["Allele_comparison"],
                "P8DS1_n_over_N_percent": row["P8DS1_n_over_N_percent"],
                "Control_group": row["Control_group"],
                "Control_n_over_N_percent": row["Control_n_over_N_percent"],
                "Odds_ratio": format_odds_ratio(row["Odds_ratio"]),
                "Fisher_P_raw": format_p(row["Fisher_P_raw"]),
                "BH_adjusted_P": format_p(row["BH_adjusted_P"]),
                "Significant_raw": row["Significant_raw"],
                "Significant_after_BH": row["Significant_after_BH"],
            })

    print("\nFisher analysis complete.")
    print(f"Manuscript-ready Table S4 written to:\n{OUTPUT_FILE}\n")

    print(
        "Segment\tAllele comparison\tP[8]DS-1 n/N (%)\t"
        "Control group\tControl n/N (%)\tRaw P\tBH P\tSignificant"
    )

    for row in results:
        print(
            f"{row['Segment']}\t"
            f"{row['Allele_comparison']}\t"
            f"{row['P8DS1_n_over_N_percent']}\t"
            f"{row['Control_group']}\t"
            f"{row['Control_n_over_N_percent']}\t"
            f"{format_p(row['Fisher_P_raw'])}\t"
            f"{format_p(row['BH_adjusted_P'])}\t"
            f"{row['Significant_after_BH']}"
        )


if __name__ == "__main__":
    main()
