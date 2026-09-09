#!/usr/bin/env python3

import csv
import math
import re
from pathlib import Path

SEGMENTS = [
    "VP1", "VP2", "VP3", "VP4", "VP6",
    "NSP1", "NSP2", "NSP3", "NSP4"
]

BASE_DIR = Path("/deac/bio/esstmanGrp/otool23/P1_Genetic distances")
OUTPUT = BASE_DIR / "Allele_identity_distance_table.csv"


def allele(label):
    m = re.match(r"^\s*([ABCD])-", str(label), flags=re.I)
    return m.group(1).upper() if m else None


def read_identity_matrix(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))

    if len(rows) < 2:
        raise ValueError(f"{path}: matrix is empty or incomplete")

    col_labels = [x.strip() for x in rows[0][1:]]
    row_labels = [r[0].strip() for r in rows[1:]]

    matrix = []

    for r in rows[1:]:
        vals = []

        for x in r[1:]:
            x = x.strip()
            vals.append(
                float(x)
                if x not in ("", "NA", "NaN", "nan")
                else math.nan
            )

        matrix.append(vals)

    if len(matrix) != len(row_labels):
        raise ValueError(f"Row-label mismatch in {path}")

    if any(len(r) != len(col_labels) for r in matrix):
        raise ValueError(f"Column mismatch in {path}")

    if len(row_labels) != len(col_labels):
        raise ValueError(f"Matrix is not square in {path}")

    if row_labels != col_labels:
        raise ValueError(
            f"Row and column labels are not identical/in the same order in {path}"
        )

    return row_labels, matrix


def mean(values):
    values = [v for v in values if not math.isnan(v)]
    return sum(values) / len(values) if values else math.nan


def fmt(value):
    if value is None or math.isnan(value):
        return ""
    return f"{value:.4f}"


def summarize_segment(segment, path):
    labels, identity = read_identity_matrix(path)

    groups = {
        g: [
            i for i, lab in enumerate(labels)
            if allele(lab) == g
        ]
        for g in "ABCD"
    }

    print(
        f"{segment}: "
        f"A={len(groups['A'])}, "
        f"B={len(groups['B'])}, "
        f"C={len(groups['C'])}, "
        f"D={len(groups['D'])}"
    )

    def identity_value(i, j):
        return identity[i][j]

    def distance_value(i, j):
        v = identity[i][j]
        return math.nan if math.isnan(v) else 100.0 - v

    def within_identity(indices):
        return [
            identity_value(indices[i], indices[j])
            for i in range(len(indices))
            for j in range(i + 1, len(indices))
        ]

    def within_distance(indices):
        return [
            distance_value(indices[i], indices[j])
            for i in range(len(indices))
            for j in range(i + 1, len(indices))
        ]

    def between_identity(indices1, indices2):
        return [
            identity_value(i, j)
            for i in indices1
            for j in indices2
        ]

    def between_distance(indices1, indices2):
        return [
            distance_value(i, j)
            for i in indices1
            for j in indices2
        ]

    row = {"Segment": segment}

    within_identity_means = {}
    within_distance_means = {}

    for g in "ABCD":
        indices = groups[g]

        within_id = mean(within_identity(indices))
        within_dist = mean(within_distance(indices))

        within_identity_means[g] = within_id
        within_distance_means[g] = within_dist

        row[f"n_{g}"] = len(indices)
        row[f"Within_{g}_identity_%"] = fmt(within_id)
        row[f"Within_{g}_distance_%"] = fmt(within_dist)

    B = groups["B"]

    for other in ("A", "C", "D"):
        G = groups[other]
        prefix = f"B_vs_{other}"

        if not B or not G:
            row[f"{prefix}_between_identity_%"] = ""
            row[f"{prefix}_between_distance_%"] = ""
            row[f"{prefix}_identity_ratio"] = ""
            row[f"{prefix}_distance_ratio"] = ""
            continue

        between_id = mean(
            between_identity(B, G)
        )

        between_dist = mean(
            between_distance(B, G)
        )

        identity_ratio = (
            within_identity_means["B"] / between_id
            if (
                not math.isnan(within_identity_means["B"])
                and not math.isnan(between_id)
                and between_id != 0
            )
            else math.nan
        )

        distance_ratio = (
            within_distance_means["B"] / between_dist
            if (
                not math.isnan(within_distance_means["B"])
                and not math.isnan(between_dist)
                and between_dist != 0
            )
            else math.nan
        )

        row[f"{prefix}_between_identity_%"] = fmt(between_id)
        row[f"{prefix}_between_distance_%"] = fmt(between_dist)
        row[f"{prefix}_identity_ratio"] = fmt(identity_ratio)
        row[f"{prefix}_distance_ratio"] = fmt(distance_ratio)

    return row


all_rows = []

for segment in SEGMENTS:
    path = BASE_DIR / f"{segment}.csv"

    if not path.exists():
        print(f"Not found: {path.name} — skipped.")
        continue

    all_rows.append(
        summarize_segment(segment, path)
    )


FIELDS = [
    "Segment",

    "n_A",
    "Within_A_identity_%",
    "Within_A_distance_%",

    "n_B",
    "Within_B_identity_%",
    "Within_B_distance_%",

    "n_C",
    "Within_C_identity_%",
    "Within_C_distance_%",

    "n_D",
    "Within_D_identity_%",
    "Within_D_distance_%",

    "B_vs_A_between_identity_%",
    "B_vs_A_between_distance_%",
    "B_vs_A_identity_ratio",
    "B_vs_A_distance_ratio",

    "B_vs_C_between_identity_%",
    "B_vs_C_between_distance_%",
    "B_vs_C_identity_ratio",
    "B_vs_C_distance_ratio",

    "B_vs_D_between_identity_%",
    "B_vs_D_between_distance_%",
    "B_vs_D_identity_ratio",
    "B_vs_D_distance_ratio",
]

with open(
    OUTPUT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=FIELDS
    )

    writer.writeheader()

    for row in all_rows:
        writer.writerow(row)


print(f"Done. Results written to: {OUTPUT}")
