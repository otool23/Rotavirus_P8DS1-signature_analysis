#!/usr/bin/env python3
"""
Paper 1: amino-acid and nucleotide signature discovery
=======================================================

Purpose
-------
Identify accepted amino-acid (AA) and nucleotide (NU) signatures associated
with the P[8]DS-1 reassortant / type-B background using callable sites only.

This streamlined version produces ONLY TWO output files:

    AA_signature_changes.csv
    NU_signature_changes.csv

No Fisher's exact tests, graphs, heatmaps, candidate-audit files,
failure-detail files, metadata files, or other outputs are generated.

Input files
-----------
VP1_AA.fasta   VP1_NU.fasta
VP2_AA.fasta   VP2_NU.fasta
VP3_AA.fasta   VP3_NU.fasta
VP4_AA.fasta   VP4_NU.fasta
VP6_AA.fasta   VP6_NU.fasta
NSP1_AA.fasta  NSP1_NU.fasta
NSP2_AA.fasta  NSP2_NU.fasta
NSP3_AA.fasta  NSP3_NU.fasta
NSP4_AA.fasta  NSP4_NU.fasta

VP7 and NSP5 are intentionally excluded.

Header-prefix groups
--------------------
Non-VP4:
    RP_   P[8]DS-1 reassortant, type-B/pink allele
    WTP_  wild-type strain carrying type-B/pink allele
    R_    type-A/red allele
    O_    type-C/orange allele
    G_    type-D/grey allele

VP4:
    LGR_  reassortant/light-green group
    LGWT_ wild-type strain carrying light-green/reassortant-like VP4
    DG_   wild-type/deep-green group

Signature criteria
------------------
Criterion 1:
    Non-VP4:
        RP  = 100% among callable sequences
        WTP = 100% among callable sequences
        R/O/G = 0% among callable sequences, if the group exists

    VP4:
        LGR  = 100%
        LGWT = 100%
        DG   = 0%

Criterion 2:
    Non-VP4:
        RP  = 100% among callable sequences
        0% < WTP < 100%
        R/O/G = 0% among callable sequences, if the group exists

    VP4:
        LGR  = 100%
        0% < LGWT < 100%
        DG   = 0%

Callable-site rule
------------------
Gaps and ambiguous symbols are excluded from site-specific denominators.

    prevalence (%) = n_state / callable_N * 100

An existing required group must meet the minimum callable fraction
(default = 0.80). A group absent entirely is treated as not applicable.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


SCRIPT_VERSION = "3.0.0"

SEGMENTS = [
    "VP1", "VP2", "VP3", "VP4", "VP6",
    "NSP1", "NSP2", "NSP3", "NSP4",
]

SEQ_TYPES = ["AA", "NU"]

NON_VP4_GROUPS = ["RP", "WTP", "R", "O", "G"]
VP4_GROUPS = ["LGR", "LGWT", "DG"]

AA_VALID = set("ACDEFGHIKLMNPQRSTVWY")
NU_VALID = set("ACGTU")


# ============================================================================
# COMMAND-LINE ARGUMENTS
# ============================================================================

def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description=(
            "Identify accepted AA and nucleotide signatures using callable "
            "sites only. Produces exactly two CSV output files."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=script_dir,
        help=(
            "Directory containing *_AA.fasta and *_NU.fasta files "
            "(default: directory containing this script)."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "P1_signature_results",
        help=(
            "Directory for the two CSV outputs "
            "(default: P1_signature_results beside the script)."
        ),
    )

    parser.add_argument(
        "--min-callable-fraction",
        type=float,
        default=0.80,
        help=(
            "Minimum callable fraction required for an existing group "
            "(default: 0.80)."
        ),
    )

    args = parser.parse_args()

    if not 0.0 < args.min_callable_fraction <= 1.0:
        parser.error("--min-callable-fraction must be >0 and <=1.")

    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()

    return args


# ============================================================================
# FASTA / GROUPING
# ============================================================================

def classify_header(header: str, segment: str) -> str | None:
    name = header.strip().upper()

    if segment == "VP4":
        for prefix in ("LGWT_", "LGR_", "DG_"):
            if name.startswith(prefix):
                return prefix[:-1]
        return None

    for prefix in ("WTP_", "RP_", "R_", "O_", "G_"):
        if name.startswith(prefix):
            return prefix[:-1]

    return None


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    header = None
    pieces: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    if header in sequences:
                        raise ValueError(
                            f"Duplicate FASTA header in {path.name}: {header}"
                        )
                    sequences[header] = "".join(pieces).replace(" ", "").upper()

                header = line[1:].strip()
                pieces = []

            else:
                if header is None:
                    raise ValueError(
                        f"{path.name}: sequence encountered before FASTA header."
                    )
                pieces.append(line)

        if header is not None:
            if header in sequences:
                raise ValueError(
                    f"Duplicate FASTA header in {path.name}: {header}"
                )
            sequences[header] = "".join(pieces).replace(" ", "").upper()

    if not sequences:
        raise ValueError(f"No sequences found in {path}")

    lengths = {len(seq) for seq in sequences.values()}

    if len(lengths) != 1:
        raise ValueError(
            f"{path.name} is not a valid multiple-sequence alignment: "
            f"sequence lengths differ ({sorted(lengths)})."
        )

    return sequences


def split_groups(
    sequences: dict[str, str],
    segment: str,
) -> dict[str, list[str]]:

    expected = VP4_GROUPS if segment == "VP4" else NON_VP4_GROUPS
    grouped = {group: [] for group in expected}

    for header in sequences:
        group = classify_header(header, segment)

        if group in grouped:
            grouped[group].append(header)

    for group in grouped:
        grouped[group].sort()

    return grouped


# ============================================================================
# COORDINATE / CALLABLE-STATE HELPERS
# ============================================================================

def ungapped_length(sequence: str) -> int:
    return sum(char != "-" for char in sequence)


def choose_coordinate_reference(
    sequences: dict[str, str],
    segment: str,
) -> tuple[str, str]:

    wt_group = "DG" if segment == "VP4" else "R"

    candidates = [
        (header, seq)
        for header, seq in sequences.items()
        if classify_header(header, segment) == wt_group
    ]

    if not candidates:
        raise ValueError(
            f"{segment}: no {wt_group}_ sequence available for coordinate numbering."
        )

    candidates.sort(
        key=lambda item: (-ungapped_length(item[1]), item[0])
    )

    return candidates[0]


def build_position_map(reference_sequence: str) -> list[int | None]:
    mapping: list[int | None] = []
    position = 0

    for char in reference_sequence:
        if char == "-":
            mapping.append(None)
        else:
            position += 1
            mapping.append(position)

    return mapping


def is_valid_state(char: str, seq_type: str) -> bool:
    if seq_type == "AA":
        return char in AA_VALID
    return char in NU_VALID


def modal_state(chars: list[str], seq_type: str) -> str | None:
    valid = [
        char
        for char in chars
        if is_valid_state(char, seq_type)
    ]

    if not valid:
        return None

    ranked = Counter(valid).most_common()

    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return None

    return ranked[0][0]


def group_chars(
    headers: list[str],
    sequences: dict[str, str],
    column: int,
) -> list[str]:

    return [
        sequences[header][column]
        for header in headers
    ]


def state_stats(
    headers: list[str],
    sequences: dict[str, str],
    column: int,
    state: str,
    seq_type: str,
) -> dict[str, float | int]:

    values = group_chars(
        headers,
        sequences,
        column,
    )

    total_n = len(values)

    callable_values = [
        value
        for value in values
        if is_valid_state(value, seq_type)
    ]

    callable_n = len(callable_values)
    missing_n = total_n - callable_n
    n_state = sum(value == state for value in callable_values)

    pct = (
        math.nan
        if callable_n == 0
        else 100.0 * n_state / callable_n
    )

    callable_fraction = (
        math.nan
        if total_n == 0
        else callable_n / total_n
    )

    return {
        "n": n_state,
        "total_N": total_n,
        "callable_N": callable_n,
        "missing_N": missing_n,
        "pct": pct,
        "callable_fraction": callable_fraction,
    }


def coverage_is_sufficient(
    stat: dict[str, float | int],
    min_callable_fraction: float,
) -> bool:

    if stat["total_N"] == 0:
        return False

    if stat["callable_N"] == 0:
        return False

    return (
        float(stat["callable_fraction"])
        >= min_callable_fraction
    )


def exclusion_group_passes(
    stat: dict[str, float | int],
    min_callable_fraction: float,
) -> bool:

    # Completely absent group = not applicable, therefore does not fail.
    if stat["total_N"] == 0:
        return True

    # Existing group must have adequate callable coverage.
    if not coverage_is_sufficient(
        stat,
        min_callable_fraction,
    ):
        return False

    # Candidate state must be absent from the exclusion group.
    return stat["pct"] == 0.0


# ============================================================================
# OUTPUT TABLE
# ============================================================================

GROUP_ROLES = [
    "Reassortant",
    "WT_B",
    "WT_A",
    "C",
    "D",
]

SIGNATURE_FIELDS = [
    "Segment",
    "Position",
    "Alignment_column",
    "Mutation",
    "WT_state",
    "Reassortant_state",
    "Criterion",
    "Coordinate_reference",
    "Min_callable_fraction",
]

for role in GROUP_ROLES:
    SIGNATURE_FIELDS.extend([
        f"{role}_group",
        f"{role}_n",
        f"{role}_total_N",
        f"{role}_callable_N",
        f"{role}_missing_N",
        f"{role}_callable_fraction",
        f"{role}_pct",
        f"{role}_n_over_callable",
    ])


def blank_stat_fields(role: str) -> dict[str, object]:
    return {
        f"{role}_group": "",
        f"{role}_n": "",
        f"{role}_total_N": "",
        f"{role}_callable_N": "",
        f"{role}_missing_N": "",
        f"{role}_callable_fraction": "",
        f"{role}_pct": "",
        f"{role}_n_over_callable": "",
    }


def add_group_stat(
    row: dict[str, object],
    role: str,
    group_name: str,
    stat: dict[str, float | int],
) -> None:

    row[f"{role}_group"] = group_name
    row[f"{role}_n"] = stat["n"]
    row[f"{role}_total_N"] = stat["total_N"]
    row[f"{role}_callable_N"] = stat["callable_N"]
    row[f"{role}_missing_N"] = stat["missing_N"]

    if stat["total_N"] == 0:
        row[f"{role}_callable_fraction"] = ""
        row[f"{role}_pct"] = ""
        row[f"{role}_n_over_callable"] = "NA"
        return

    row[f"{role}_callable_fraction"] = stat["callable_fraction"]
    row[f"{role}_pct"] = stat["pct"]

    if stat["callable_N"] == 0:
        row[f"{role}_n_over_callable"] = f"{stat['n']}/0"
    else:
        row[f"{role}_n_over_callable"] = (
            f"{stat['n']}/{stat['callable_N']}"
        )


def build_site_row(
    *,
    segment: str,
    position: int,
    alignment_column: int,
    mutation: str,
    wt_state: str,
    reassortant_state: str,
    criterion: str,
    coordinate_reference: str,
    min_callable_fraction: float,
    stats: dict[str, dict[str, float | int]],
) -> dict[str, object]:

    row: dict[str, object] = {
        "Segment": segment,
        "Position": position,
        "Alignment_column": alignment_column,
        "Mutation": mutation,
        "WT_state": wt_state,
        "Reassortant_state": reassortant_state,
        "Criterion": criterion,
        "Coordinate_reference": coordinate_reference,
        "Min_callable_fraction": min_callable_fraction,
    }

    for role in GROUP_ROLES:
        row.update(
            blank_stat_fields(role)
        )

    if segment == "VP4":
        mapping = {
            "Reassortant": "LGR",
            "WT_B": "LGWT",
            "WT_A": "DG",
        }

    else:
        mapping = {
            "Reassortant": "RP",
            "WT_B": "WTP",
            "WT_A": "R",
            "C": "O",
            "D": "G",
        }

    for role, group_name in mapping.items():
        add_group_stat(
            row,
            role,
            group_name,
            stats[group_name],
        )

    return row


# ============================================================================
# CORE SIGNATURE ANALYSIS
# ============================================================================

def analyze_alignment(
    segment: str,
    seq_type: str,
    fasta_path: Path,
    min_callable_fraction: float,
) -> list[dict[str, object]]:

    sequences = read_fasta(
        fasta_path
    )

    grouped = split_groups(
        sequences,
        segment,
    )

    expected_groups = (
        VP4_GROUPS
        if segment == "VP4"
        else NON_VP4_GROUPS
    )

    if segment == "VP4":
        reassortant_group = "LGR"
        wt_b_group = "LGWT"
        wt_reference_group = "DG"
        exclusion_groups = ["DG"]

    else:
        reassortant_group = "RP"
        wt_b_group = "WTP"
        wt_reference_group = "R"
        exclusion_groups = ["R", "O", "G"]

    missing_required_groups = [
        group
        for group in (
            reassortant_group,
            wt_b_group,
            wt_reference_group,
        )
        if len(grouped[group]) == 0
    ]

    if missing_required_groups:
        print(
            f"{segment} {seq_type}: skipped; missing required group(s): "
            + ", ".join(missing_required_groups)
        )
        return []

    reference_header, reference_sequence = (
        choose_coordinate_reference(
            sequences,
            segment,
        )
    )

    position_map = build_position_map(
        reference_sequence
    )

    accepted: list[dict[str, object]] = []

    for col in range(
        len(reference_sequence)
    ):

        position = position_map[col]

        if position is None:
            continue

        # Candidate reassortant state.
        reassortant_chars = group_chars(
            grouped[reassortant_group],
            sequences,
            col,
        )

        reassortant_state = modal_state(
            reassortant_chars,
            seq_type,
        )

        if reassortant_state is None:
            continue

        reass_stat = state_stats(
            grouped[reassortant_group],
            sequences,
            col,
            reassortant_state,
            seq_type,
        )

        # Reassortant group must have sufficient callable data.
        if not coverage_is_sufficient(
            reass_stat,
            min_callable_fraction,
        ):
            continue

        # Reassortant state must be fixed among callable reassortants.
        if reass_stat["pct"] != 100.0:
            continue

        # Determine WT state from R_ or DG_.
        wt_chars = group_chars(
            grouped[wt_reference_group],
            sequences,
            col,
        )

        wt_state = modal_state(
            wt_chars,
            seq_type,
        )

        if wt_state is None:
            continue

        if wt_state == reassortant_state:
            continue

        mutation = (
            f"{wt_state}{position}{reassortant_state}"
        )

        stats = {
            group: state_stats(
                grouped[group],
                sequences,
                col,
                reassortant_state,
                seq_type,
            )
            for group in expected_groups
        }

        wt_b_stat = stats[wt_b_group]

        # WT-B group must also meet callable threshold.
        wt_b_coverage_ok = coverage_is_sufficient(
            wt_b_stat,
            min_callable_fraction,
        )

        if not wt_b_coverage_ok:
            continue

        # All existing exclusion groups must have adequate callable data
        # and 0% of the reassortant-associated state.
        exclusions_pass = all(
            exclusion_group_passes(
                stats[group],
                min_callable_fraction,
            )
            for group in exclusion_groups
        )

        if not exclusions_pass:
            continue

        criterion1 = (
            wt_b_stat["pct"] == 100.0
        )

        criterion2 = (
            not math.isnan(
                float(wt_b_stat["pct"])
            )
            and 0.0
            < float(wt_b_stat["pct"])
            < 100.0
        )

        if criterion1:
            criterion = (
                "Criterion_1_Strict_B_signature"
            )

        elif criterion2:
            criterion = (
                "Criterion_2_Reassortant_enriched_B_signature"
            )

        else:
            continue

        accepted.append(
            build_site_row(
                segment=segment,
                position=position,
                alignment_column=col + 1,
                mutation=mutation,
                wt_state=wt_state,
                reassortant_state=reassortant_state,
                criterion=criterion,
                coordinate_reference=reference_header,
                min_callable_fraction=min_callable_fraction,
                stats=stats,
            )
        )

    print(
        f"{segment} {seq_type}: "
        f"{len(accepted)} accepted signatures"
    )

    return accepted


# ============================================================================
# CSV WRITING
# ============================================================================

def format_output_value(
    value: object,
) -> object:

    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6f}"

    return value


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=SIGNATURE_FIELDS,
        )

        writer.writeheader()

        for row in rows:
            writer.writerow({
                key: format_output_value(
                    row.get(key, "")
                )
                for key in SIGNATURE_FIELDS
            })


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    args = parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    aa_signatures: list[dict[str, object]] = []
    nu_signatures: list[dict[str, object]] = []

    for segment in SEGMENTS:
        for seq_type in SEQ_TYPES:

            fasta_path = (
                args.input_dir
                / f"{segment}_{seq_type}.fasta"
            )

            if not fasta_path.exists():
                print(
                    f"Not found: {fasta_path.name} — skipped."
                )
                continue

            signatures = analyze_alignment(
                segment=segment,
                seq_type=seq_type,
                fasta_path=fasta_path,
                min_callable_fraction=(
                    args.min_callable_fraction
                ),
            )

            if seq_type == "AA":
                aa_signatures.extend(
                    signatures
                )
            else:
                nu_signatures.extend(
                    signatures
                )

    segment_order = {
        segment: i
        for i, segment in enumerate(SEGMENTS)
    }

    def sort_key(
        row: dict[str, object],
    ):
        return (
            segment_order[
                str(row["Segment"])
            ],
            int(row["Position"]),
            str(row["Mutation"]),
        )

    aa_signatures.sort(
        key=sort_key
    )

    nu_signatures.sort(
        key=sort_key
    )

    aa_output = (
        args.output_dir
        / "AA_signature_changes.csv"
    )

    nu_output = (
        args.output_dir
        / "NU_signature_changes.csv"
    )

    write_csv(
        aa_output,
        aa_signatures,
    )

    write_csv(
        nu_output,
        nu_signatures,
    )

    print("\n" + "=" * 72)
    print("ANALYSIS COMPLETE")
    print("=" * 72)
    print(
        f"Accepted AA signatures: "
        f"{len(aa_signatures)}"
    )
    print(
        f"Accepted NU signatures: "
        f"{len(nu_signatures)}"
    )
    print(
        f"AA CSV: {aa_output}"
    )
    print(
        f"NU CSV: {nu_output}"
    )
    print(
        "\nOnly the two signature CSV files were generated."
    )


if __name__ == "__main__":
    main()
