#!/usr/bin/env python3
"""
Build one concatenated nucleotide MSA for reassortant Rotavirus A strains,
with automatic NCBI fallback for missing local accessions.

Workflow
--------
1. Read a strain-to-segment accession table.
2. For each required segment:
   a. Search the local segment FASTA for the listed accession.
   b. If the accession is missing locally, download that exact accession
      from NCBI nucleotide.
   c. Verify that the downloaded record is consistent with the expected
      Rotavirus A segment/gene.
   d. Save downloaded fallback sequences separately without modifying the
      original local FASTA files.
3. Require every retained strain to have all required segments.
4. Align each segment separately with MAFFT.
5. Reorder every segment MSA so strains occur in exactly the same order.
6. Concatenate the aligned segments strain-by-strain.
7. Write ONE final concatenated FASTA MSA, where each FASTA record is one strain.
8. Calculate pairwise percent identity and genetic distance from the concatenated MSA.
9. Write square pairwise identity and distance matrices.
10. Summarize mean within-clade and between-clade identity and genetic distance.

Required segments
-----------------
VP4, VP6, VP1, VP2, VP3, NSP1, NSP2, NSP3, NSP4

Expected local nucleotide FASTA files
-------------------------------------
VP4_NU.fasta
VP6_NU.fasta
VP1_NU.fasta
VP2_NU.fasta
VP3_NU.fasta
NSP1_NU.fasta
NSP2_NU.fasta
NSP3_NU.fasta
NSP4_NU.fasta

Recognized source FASTA prefixes
--------------------------------
Non-VP4:
    RP_
    WTP_
    R_
    O_
    G_

VP4:
    LGR_
    LGWT_
    DG_

The accession is extracted from the FASTA header after the prefix.

Accession table
---------------
Required CSV columns:

    Clade,VP4,VP6,VP1,VP2,VP3,NSP1,NSP2,NSP3,NSP4

Optional:
    Strain

NCBI fallback
-------------
If an accession is absent from the corresponding local FASTA, the script:
- downloads the exact nucleotide record from NCBI;
- verifies the record against the expected segment/gene;
- saves it under:
      <output-dir>/00_downloaded_missing_sequences/<SEGMENT>_NCBI_missing.fasta
- uses it together with the local sequences for MAFFT and concatenation.

Original local FASTA files are never modified.

MAFFT
-----
MAFFT must be installed and available in PATH.

Example:
    module load apps/mafft/7.525

Run
---
module load apps/mafft/7.525

/deac/bio/esstmanGrp/otool23/ml_env/bin/python \
"/deac/bio/esstmanGrp/otool23/clades/Constellation_clades_definition.py" \
--accession-table "/deac/bio/esstmanGrp/otool23/clades/Clades_accessions.csv" \
--input-dir "/deac/bio/esstmanGrp/otool23/clades" \
--output-dir "/deac/bio/esstmanGrp/otool23/clades/concatenated_MSA" \
--ncbi-email otool23@wfu.edu

Main outputs
------------
reassortant_concatenated_MSA.fasta
clades_pairwise_identity_matrix.csv
clades_pairwise_genetic_distance_matrix.csv
clade_within_between_identity_distance_summary.csv

Each FASTA record is ONE reassortant strain, with aligned segments concatenated
in this fixed order:

    VP4 | VP6 | VP1 | VP2 | VP3 | NSP1 | NSP2 | NSP3 | NSP4
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import sys
import time
from collections import OrderedDict, defaultdict
from pathlib import Path


SCRIPT_VERSION = "3.0.0"

SEGMENTS = [
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

NON_VP4_PREFIXES = ("RP_", "WTP_", "R_", "O_", "G_")
VP4_PREFIXES = ("LGR_", "LGWT_", "DG_")

ACCESSION_RE = re.compile(r"([A-Z]{1,4}_?\d+(?:\.\d+)?)", re.I)

# Terms used to check that an NCBI record matches the expected RVA segment/gene.
SEGMENT_TERMS = {
    "VP4": ["VP4", "SEGMENT 4", "GENOME SEGMENT 4"],
    "VP6": ["VP6", "SEGMENT 6", "GENOME SEGMENT 6"],
    "VP1": ["VP1", "SEGMENT 1", "GENOME SEGMENT 1"],
    "VP2": ["VP2", "SEGMENT 2", "GENOME SEGMENT 2"],
    "VP3": ["VP3", "SEGMENT 3", "GENOME SEGMENT 3"],
    "NSP1": ["NSP1", "SEGMENT 5", "GENOME SEGMENT 5"],
    "NSP2": ["NSP2", "SEGMENT 8", "GENOME SEGMENT 8"],
    "NSP3": ["NSP3", "SEGMENT 7", "GENOME SEGMENT 7"],
    "NSP4": ["NSP4", "SEGMENT 10", "GENOME SEGMENT 10"],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a concatenated RVA reassortant MSA with automatic NCBI "
            "fallback for accessions missing from local segment FASTA files."
        )
    )

    parser.add_argument(
        "--accession-table",
        type=Path,
        required=True,
        help="CSV containing Clade and segment accession columns.",
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("."),
        help="Directory containing *_NU.fasta files.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("concatenated_MSA"),
        help="Output directory.",
    )

    parser.add_argument(
        "--mafft",
        default="mafft",
        help="MAFFT executable name/path (default: mafft).",
    )

    parser.add_argument(
        "--mafft-threads",
        default="-1",
        help="MAFFT thread count; -1 lets MAFFT choose automatically.",
    )

    parser.add_argument(
        "--ncbi-email",
        required=True,
        help="Email required by NCBI Entrez.",
    )

    parser.add_argument(
        "--ncbi-api-key",
        default=None,
        help="Optional NCBI API key.",
    )

    parser.add_argument(
        "--no-segment-check",
        action="store_true",
        help=(
            "Disable NCBI segment/gene verification. Not recommended unless "
            "a valid record cannot be recognized by the automatic checker."
        ),
    )

    return parser.parse_args()


# ============================================================================
# FASTA UTILITIES
# ============================================================================

def read_fasta(path: Path):
    records = OrderedDict()
    header = None
    seq_parts = []

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    if header in records:
                        raise ValueError(
                            f"Duplicate FASTA header in {path.name}: {header}"
                        )
                    records[header] = "".join(seq_parts).replace(" ", "").upper()

                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(
                        f"{path.name}: sequence encountered before FASTA header."
                    )
                seq_parts.append(line)

    if header is not None:
        if header in records:
            raise ValueError(
                f"Duplicate FASTA header in {path.name}: {header}"
            )
        records[header] = "".join(seq_parts).replace(" ", "").upper()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def write_fasta(records, path: Path, line_width=80):
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records.items():
            handle.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                handle.write(seq[i:i + line_width] + "\n")


def append_fasta_record(path: Path, header: str, sequence: str, line_width=80):
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            handle.write(sequence[i:i + line_width] + "\n")


def accession_from_header(header: str, segment: str):
    prefixes = VP4_PREFIXES if segment == "VP4" else NON_VP4_PREFIXES

    for prefix in prefixes:
        if header.startswith(prefix):
            remainder = header[len(prefix):]
            token = re.split(r"[\s|;/]", remainder, maxsplit=1)[0].strip()

            match = ACCESSION_RE.search(token)
            if match:
                return match.group(1).upper().split(".")[0]

    match = ACCESSION_RE.search(header)
    if match:
        return match.group(1).upper().split(".")[0]

    return None


def build_accession_index(records, segment: str):
    index = {}

    for header, sequence in records.items():
        accession = accession_from_header(header, segment)

        if accession is None:
            continue

        accession = accession.upper().split(".")[0]

        if accession in index:
            raise ValueError(
                f"{segment}: accession {accession} appears more than once "
                f"in the FASTA file."
            )

        index[accession] = (header, sequence)

    return index


# ============================================================================
# ACCESSION TABLE
# ============================================================================

def read_accession_table(path: Path):
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError("Accession table has no header row.")

        required = ["Clade"] + SEGMENTS
        missing_columns = [
            col for col in required
            if col not in reader.fieldnames
        ]

        if missing_columns:
            raise ValueError(
                "Accession table is missing required column(s): "
                + ", ".join(missing_columns)
            )

        rows = []

        for row_number, row in enumerate(reader, start=2):
            cleaned = {
                key: (value.strip() if value else "")
                for key, value in row.items()
            }

            if not cleaned["Clade"]:
                raise ValueError(
                    f"Row {row_number}: Clade is missing."
                )

            missing_accessions = [
                segment for segment in SEGMENTS
                if not cleaned[segment]
            ]

            if missing_accessions:
                raise ValueError(
                    f"Row {row_number}: missing accession(s) for "
                    + ", ".join(missing_accessions)
                )

            for segment in SEGMENTS:
                cleaned[segment] = cleaned[segment].upper().split(".")[0]

            cleaned["_row_number"] = row_number
            rows.append(cleaned)

    if not rows:
        raise ValueError("Accession table contains no strains.")

    return rows


# ============================================================================
# NCBI HELPERS
# ============================================================================

def setup_entrez(email, api_key=None):
    try:
        from Bio import Entrez
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required. Install with: pip install biopython"
        ) from exc

    Entrez.email = email

    if api_key:
        Entrez.api_key = api_key

    return Entrez


def fetch_ncbi_record(accession, email, api_key=None):
    """
    Fetch exact nucleotide accession as a GenBank record.
    """
    try:
        from Bio import Entrez, SeqIO
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required. Install with: pip install biopython"
        ) from exc

    Entrez.email = email

    if api_key:
        Entrez.api_key = api_key

    with Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text",
    ) as handle:
        record = SeqIO.read(handle, "genbank")

    return record


def record_text_for_segment_check(record):
    parts = [
        str(record.id),
        str(record.name),
        str(record.description),
    ]

    for feature in record.features:
        if feature.type in ("source", "gene", "CDS"):
            for key, values in feature.qualifiers.items():
                parts.append(str(key))
                parts.extend(str(v) for v in values)

    return " ".join(parts).upper()


def verify_ncbi_segment(record, expected_segment: str):
    """
    Returns (True, evidence) if the NCBI record contains expected segment/gene
    evidence. Otherwise returns (False, diagnostic text).
    """
    text = record_text_for_segment_check(record)

    # Confirm rotavirus context when recognizable.
    rotavirus_ok = (
        "ROTAVIRUS" in text
        or "ROTAVIRUS A" in text
    )

    expected_terms = SEGMENT_TERMS[expected_segment]
    matched_terms = [
        term
        for term in expected_terms
        if term.upper() in text
    ]

    if rotavirus_ok and matched_terms:
        return True, "; ".join(matched_terms)

    return False, (
        f"Expected {expected_segment}; no convincing matching annotation "
        f"found. Record description: {record.description}"
    )


def get_ncbi_strain_name(record):
    for feature in record.features:
        if feature.type == "source":
            strain = feature.qualifiers.get("strain")
            if strain:
                return strain[0].strip()

            isolate = feature.qualifiers.get("isolate")
            if isolate:
                return isolate[0].strip()

    description = (record.description or "").strip()
    return description if description else None


def sanitize_header_text(text: str):
    text = text.strip()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.\-\[\]]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


# ============================================================================
# MAFFT
# ============================================================================

def run_mafft(input_fasta: Path, output_fasta: Path, mafft_exe: str, threads: str):
    command = [
        mafft_exe,
        "--auto",
        "--thread",
        str(threads),
        str(input_fasta),
    ]

    print("Running MAFFT:", " ".join(command))

    with output_fasta.open("w", encoding="utf-8") as out_handle:
        result = subprocess.run(
            command,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            text=True,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"MAFFT failed for {input_fasta.name}:\n{result.stderr}"
        )


# ============================================================================
# CONCATENATED-MSA IDENTITY / DISTANCE UTILITIES
# ============================================================================

CALLABLE_NT = set("ACGTU")


def pairwise_identity_percent(seq1: str, seq2: str):
    """
    Calculate pairwise nucleotide identity (%) using only sites at which BOTH
    sequences contain an unambiguous nucleotide (A/C/G/T/U).

    Gaps and ambiguous bases are excluded from the denominator. U is treated
    as equivalent to T.

    Returns:
        (identity_percent, callable_sites)
    """
    if len(seq1) != len(seq2):
        raise ValueError("Pairwise comparison requires equal-length aligned sequences.")

    matches = 0
    callable_sites = 0

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a not in CALLABLE_NT or b not in CALLABLE_NT:
            continue

        a_norm = "T" if a == "U" else a
        b_norm = "T" if b == "U" else b

        callable_sites += 1
        if a_norm == b_norm:
            matches += 1

    if callable_sites == 0:
        return math.nan, 0

    return 100.0 * matches / callable_sites, callable_sites


def build_pairwise_matrices(concatenated):
    """
    Build symmetric pairwise identity and genetic-distance matrices.

    Genetic distance (%) = 100 - percent identity.
    """
    headers = list(concatenated.keys())
    identity_matrix = {
        header: {other: math.nan for other in headers}
        for header in headers
    }
    distance_matrix = {
        header: {other: math.nan for other in headers}
        for header in headers
    }
    callable_matrix = {
        header: {other: 0 for other in headers}
        for header in headers
    }

    for i, header1 in enumerate(headers):
        for j in range(i + 1, len(headers)):
            header2 = headers[j]

            identity_pct, callable_sites = pairwise_identity_percent(
                concatenated[header1],
                concatenated[header2],
            )

            distance_pct = (
                math.nan
                if math.isnan(identity_pct)
                else 100.0 - identity_pct
            )

            identity_matrix[header1][header2] = identity_pct
            identity_matrix[header2][header1] = identity_pct

            distance_matrix[header1][header2] = distance_pct
            distance_matrix[header2][header1] = distance_pct

            callable_matrix[header1][header2] = callable_sites
            callable_matrix[header2][header1] = callable_sites

    return headers, identity_matrix, distance_matrix, callable_matrix


def write_square_matrix(path: Path, headers, matrix, decimals=4):
    """Write a square symmetric CSV matrix with a blank diagonal."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + headers)

        for row_header in headers:
            row = [row_header]

            for col_header in headers:
                if row_header == col_header:
                    row.append("")
                    continue

                value = matrix[row_header][col_header]
                row.append(
                    "" if math.isnan(value) else f"{value:.{decimals}f}"
                )

            writer.writerow(row)


def mean_nonmissing(values):
    values = [
        value for value in values
        if value is not None and not math.isnan(value)
    ]
    return sum(values) / len(values) if values else math.nan


def summarize_clades(
    headers,
    identity_matrix,
    distance_matrix,
    clade_by_header,
    clade_order,
):
    """
    Produce one table containing both within-clade and between-clade summaries.

    One row is produced for each within-clade comparison and for every unique
    pair of clades present in the accession table.
    """
    grouped = OrderedDict(
        (clade, [])
        for clade in clade_order
    )

    for header in headers:
        clade = clade_by_header[header]
        grouped.setdefault(clade, []).append(header)

    rows = []

    # Within-clade summaries.
    for clade in clade_order:
        members = grouped.get(clade, [])
        identities = []
        distances = []

        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                h1 = members[i]
                h2 = members[j]
                identities.append(identity_matrix[h1][h2])
                distances.append(distance_matrix[h1][h2])

        rows.append({
            "Comparison_type": "Within",
            "Clade_1": clade,
            "Clade_2": clade,
            "n_Clade_1": len(members),
            "n_Clade_2": len(members),
            "Pair_count": len(identities),
            "Mean_identity_%": mean_nonmissing(identities),
            "Mean_genetic_distance_%": mean_nonmissing(distances),
        })

    # Between-clade summaries.
    for i in range(len(clade_order)):
        clade1 = clade_order[i]
        members1 = grouped.get(clade1, [])

        for j in range(i + 1, len(clade_order)):
            clade2 = clade_order[j]
            members2 = grouped.get(clade2, [])

            identities = [
                identity_matrix[h1][h2]
                for h1 in members1
                for h2 in members2
            ]
            distances = [
                distance_matrix[h1][h2]
                for h1 in members1
                for h2 in members2
            ]

            rows.append({
                "Comparison_type": "Between",
                "Clade_1": clade1,
                "Clade_2": clade2,
                "n_Clade_1": len(members1),
                "n_Clade_2": len(members2),
                "Pair_count": len(identities),
                "Mean_identity_%": mean_nonmissing(identities),
                "Mean_genetic_distance_%": mean_nonmissing(distances),
            })

    return rows


def write_clade_summary(path: Path, rows):
    fields = [
        "Comparison_type",
        "Clade_1",
        "Clade_2",
        "n_Clade_1",
        "n_Clade_2",
        "Pair_count",
        "Mean_identity_%",
        "Mean_genetic_distance_%",
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            formatted = dict(row)

            for field in ("Mean_identity_%", "Mean_genetic_distance_%"):
                value = formatted[field]
                formatted[field] = (
                    "" if math.isnan(value) else f"{value:.4f}"
                )

            writer.writerow(formatted)


# ============================================================================
# MAIN
# ============================================================================

def main():
    args = parse_args()

    args.accession_table = args.accession_table.resolve()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    downloaded_dir = args.output_dir / "00_downloaded_missing_sequences"
    extracted_dir = args.output_dir / "01_extracted_segments"
    aligned_dir = args.output_dir / "02_aligned_segments"

    downloaded_dir.mkdir(exist_ok=True)
    extracted_dir.mkdir(exist_ok=True)
    aligned_dir.mkdir(exist_ok=True)

    if shutil.which(args.mafft) is None:
        raise RuntimeError(
            f"MAFFT executable '{args.mafft}' was not found in PATH."
        )

    setup_entrez(args.ncbi_email, args.ncbi_api_key)

    strain_rows = read_accession_table(args.accession_table)
    print(f"Strains in accession table: {len(strain_rows)}")

    # Cache all downloaded NCBI records so each accession is fetched only once.
    ncbi_record_cache = {}

    # ------------------------------------------------------------------------
    # Resolve final strain names using VP4 accession/record.
    # ------------------------------------------------------------------------

    final_names = []

    for row in strain_rows:
        vp4_accession = row["VP4"]

        strain_name = row.get("Strain", "").strip()

        if not strain_name:
            try:
                if vp4_accession not in ncbi_record_cache:
                    ncbi_record_cache[vp4_accession] = fetch_ncbi_record(
                        vp4_accession,
                        args.ncbi_email,
                        args.ncbi_api_key,
                    )
                    if not args.ncbi_api_key:
                        time.sleep(0.34)

                strain_name = get_ncbi_strain_name(
                    ncbi_record_cache[vp4_accession]
                )

            except Exception as exc:
                print(
                    f"WARNING: could not retrieve NCBI strain name for "
                    f"{vp4_accession}: {exc}",
                    file=sys.stderr,
                )

        if not strain_name:
            strain_name = vp4_accession

        clade = sanitize_header_text(row["Clade"])
        strain_name = sanitize_header_text(strain_name)

        final_header = f"Clade_{clade}|{strain_name}"

        if final_header in final_names:
            final_header += f"|VP4_{vp4_accession}"

        final_names.append(final_header)

    # ------------------------------------------------------------------------
    # Process segments.
    # ------------------------------------------------------------------------

    source_report_rows = []
    download_report_rows = []
    aligned_by_segment = {}
    download_counts = defaultdict(int)

    for segment in SEGMENTS:

        fasta_path = args.input_dir / f"{segment}_NU.fasta"

        if not fasta_path.exists():
            raise FileNotFoundError(
                f"Required local FASTA file not found: {fasta_path}"
            )

        local_records = read_fasta(fasta_path)
        local_index = build_accession_index(local_records, segment)

        downloaded_fasta_path = (
            downloaded_dir / f"{segment}_NCBI_missing.fasta"
        )

        # Start fresh on each run so files reflect this run exactly.
        if downloaded_fasta_path.exists():
            downloaded_fasta_path.unlink()

        extracted = OrderedDict()

        for row, final_header in zip(strain_rows, final_names):
            accession = row[segment]

            source_type = "local"
            source_header = None
            sequence = None
            verification = "Local FASTA"

            # ------------------------------------------------------------
            # 1. Try local FASTA
            # ------------------------------------------------------------
            if accession in local_index:
                source_header, sequence = local_index[accession]

            # ------------------------------------------------------------
            # 2. Automatic NCBI fallback
            # ------------------------------------------------------------
            else:
                print(
                    f"{segment}: {accession} not found locally. "
                    "Downloading from NCBI..."
                )

                try:
                    if accession not in ncbi_record_cache:
                        ncbi_record_cache[accession] = fetch_ncbi_record(
                            accession,
                            args.ncbi_email,
                            args.ncbi_api_key,
                        )

                        if not args.ncbi_api_key:
                            time.sleep(0.34)

                    record = ncbi_record_cache[accession]

                except Exception as exc:
                    raise RuntimeError(
                        f"{segment}: failed to download accession "
                        f"{accession} from NCBI: {exc}"
                    ) from exc

                # Confirm NCBI returned the requested accession.
                returned_accession = (
                    str(record.id).upper().split(".")[0]
                )

                if returned_accession != accession:
                    raise RuntimeError(
                        f"{segment}: requested {accession}, but NCBI returned "
                        f"{record.id}."
                    )

                if not args.no_segment_check:
                    ok, verification = verify_ncbi_segment(
                        record,
                        segment,
                    )

                    if not ok:
                        raise RuntimeError(
                            f"{segment}: downloaded NCBI accession {accession} "
                            f"failed segment verification.\n{verification}\n"
                            "If you have manually confirmed the record is correct, "
                            "rerun with --no-segment-check."
                        )
                else:
                    verification = "Segment verification disabled"

                sequence = str(record.seq).upper().replace(" ", "")

                if not sequence:
                    raise RuntimeError(
                        f"{segment}: NCBI accession {accession} returned "
                        "an empty sequence."
                    )

                source_type = "NCBI_fallback"
                source_header = (
                    f"NCBI_{accession}|{sanitize_header_text(record.description)}"
                )

                append_fasta_record(
                    downloaded_fasta_path,
                    source_header,
                    sequence,
                )

                download_counts[segment] += 1

                download_report_rows.append({
                    "Segment": segment,
                    "Accession": accession,
                    "NCBI_record_id": record.id,
                    "NCBI_description": record.description,
                    "Verification": verification,
                    "Sequence_length": len(sequence),
                    "Saved_to": downloaded_fasta_path.name,
                })

            extracted[final_header] = sequence

            source_report_rows.append({
                "Final_strain_header": final_header,
                "Clade": row["Clade"],
                "Segment": segment,
                "Accession": accession,
                "Source_type": source_type,
                "Source_FASTA_or_NCBI": (
                    fasta_path.name
                    if source_type == "local"
                    else "NCBI nucleotide"
                ),
                "Source_header": source_header,
                "Verification": verification,
                "Ungapped_length": len(sequence.replace("-", "")),
            })

        if len(extracted) != len(strain_rows):
            raise RuntimeError(
                f"{segment}: extracted sequence count differs from strain count."
            )

        # ------------------------------------------------------------
        # Write exact extracted strain set used for this segment.
        # ------------------------------------------------------------
        extracted_path = extracted_dir / f"{segment}_reassortants.fasta"
        write_fasta(extracted, extracted_path)

        # ------------------------------------------------------------
        # Align this segment with MAFFT.
        # ------------------------------------------------------------
        aligned_path = aligned_dir / f"{segment}_reassortants_MSA.fasta"

        run_mafft(
            extracted_path,
            aligned_path,
            args.mafft,
            args.mafft_threads,
        )

        aligned = read_fasta(aligned_path)

        expected_headers = set(final_names)
        observed_headers = set(aligned.keys())

        if observed_headers != expected_headers:
            missing = sorted(expected_headers - observed_headers)
            extra = sorted(observed_headers - expected_headers)

            raise RuntimeError(
                f"{segment}: MAFFT output headers do not match expected strains.\n"
                f"Missing: {missing}\n"
                f"Extra: {extra}"
            )

        lengths = {len(seq) for seq in aligned.values()}

        if len(lengths) != 1:
            raise RuntimeError(
                f"{segment}: aligned sequences do not have one common length."
            )

        # Force identical strain order for all segment MSAs.
        aligned_by_segment[segment] = OrderedDict(
            (header, aligned[header])
            for header in final_names
        )

        write_fasta(
            aligned_by_segment[segment],
            aligned_path,
        )

    # ------------------------------------------------------------------------
    # Concatenate aligned segments.
    # ------------------------------------------------------------------------

    segment_lengths = {}

    for segment in SEGMENTS:
        first_sequence = next(
            iter(aligned_by_segment[segment].values())
        )
        segment_lengths[segment] = len(first_sequence)

    concatenated = OrderedDict()

    for final_header in final_names:
        concatenated[final_header] = "".join(
            aligned_by_segment[segment][final_header]
            for segment in SEGMENTS
        )

    final_lengths = {
        len(seq)
        for seq in concatenated.values()
    }

    if len(final_lengths) != 1:
        raise RuntimeError(
            "Final concatenated sequences do not have one common alignment length."
        )

    final_output = (
        args.output_dir
        / "reassortant_concatenated_MSA.fasta"
    )

    write_fasta(
        concatenated,
        final_output,
    )

    # ------------------------------------------------------------------------
    # Pairwise whole-concatenated-genome identity and genetic distance.
    # ------------------------------------------------------------------------

    clade_by_header = {
        final_header: row["Clade"]
        for row, final_header in zip(strain_rows, final_names)
    }

    # Preserve clade order as it first appears in the accession table.
    clade_order = list(OrderedDict.fromkeys(
        row["Clade"] for row in strain_rows
    ))

    (
        matrix_headers,
        identity_matrix,
        distance_matrix,
        callable_matrix,
    ) = build_pairwise_matrices(concatenated)

    identity_matrix_path = (
        args.output_dir
        / "clades_pairwise_identity_matrix.csv"
    )

    distance_matrix_path = (
        args.output_dir
        / "clades_pairwise_genetic_distance_matrix.csv"
    )

    write_square_matrix(
        identity_matrix_path,
        matrix_headers,
        identity_matrix,
        decimals=4,
    )

    write_square_matrix(
        distance_matrix_path,
        matrix_headers,
        distance_matrix,
        decimals=4,
    )

    clade_summary_rows = summarize_clades(
        matrix_headers,
        identity_matrix,
        distance_matrix,
        clade_by_header,
        clade_order,
    )

    clade_summary_path = (
        args.output_dir
        / "clade_within_between_identity_distance_summary.csv"
    )

    write_clade_summary(
        clade_summary_path,
        clade_summary_rows,
    )

    # ------------------------------------------------------------------------
    # Partition coordinates.
    # ------------------------------------------------------------------------

    partition_path = (
        args.output_dir
        / "concatenated_MSA_partitions.csv"
    )

    start = 1

    with partition_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.writer(handle)

        writer.writerow([
            "Segment",
            "Start",
            "End",
            "Aligned_length",
        ])

        for segment in SEGMENTS:
            length = segment_lengths[segment]
            end = start + length - 1

            writer.writerow([
                segment,
                start,
                end,
                length,
            ])

            start = end + 1

    # ------------------------------------------------------------------------
    # Source report.
    # ------------------------------------------------------------------------

    source_report = (
        args.output_dir
        / "concatenated_MSA_sequence_sources.csv"
    )

    with source_report.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        fieldnames = [
            "Final_strain_header",
            "Clade",
            "Segment",
            "Accession",
            "Source_type",
            "Source_FASTA_or_NCBI",
            "Source_header",
            "Verification",
            "Ungapped_length",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(source_report_rows)

    # ------------------------------------------------------------------------
    # NCBI download report.
    # ------------------------------------------------------------------------

    download_report = (
        args.output_dir
        / "NCBI_fallback_downloads.csv"
    )

    with download_report.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        fieldnames = [
            "Segment",
            "Accession",
            "NCBI_record_id",
            "NCBI_description",
            "Verification",
            "Sequence_length",
            "Saved_to",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(download_report_rows)

    # ------------------------------------------------------------------------
    # Summary.
    # ------------------------------------------------------------------------

    summary_path = (
        args.output_dir
        / "concatenated_MSA_summary.txt"
    )

    final_length = next(iter(final_lengths))
    total_downloads = sum(download_counts.values())

    with summary_path.open(
        "w",
        encoding="utf-8",
    ) as handle:

        handle.write("Paper 1 reassortant concatenated MSA\n")
        handle.write("====================================\n")
        handle.write(f"Script version: {SCRIPT_VERSION}\n")
        handle.write(f"Number of strains: {len(concatenated)}\n")
        handle.write(f"Number of segments: {len(SEGMENTS)}\n")
        handle.write(
            "Segment order: "
            + " | ".join(SEGMENTS)
            + "\n"
        )

        handle.write(
            f"Total NCBI fallback sequences downloaded: "
            f"{total_downloads}\n"
        )

        for segment in SEGMENTS:
            handle.write(
                f"{segment}: aligned length={segment_lengths[segment]}, "
                f"NCBI fallback downloads={download_counts[segment]}\n"
            )

        handle.write(
            f"Final concatenated alignment length: "
            f"{final_length}\n"
        )

        handle.write(
            f"Final MSA: {final_output.name}\n"
        )
        handle.write(
            f"Pairwise identity matrix: {identity_matrix_path.name}\n"
        )
        handle.write(
            f"Pairwise genetic-distance matrix: {distance_matrix_path.name}\n"
        )
        handle.write(
            f"Within/between clade summary: {clade_summary_path.name}\n"
        )

    print("\n" + "=" * 72)
    print("CONCATENATED MSA COMPLETE")
    print("=" * 72)

    print(
        f"Strains retained: {len(concatenated)}"
    )

    print(
        "Segment order: "
        + " | ".join(SEGMENTS)
    )

    print(
        f"NCBI fallback sequences downloaded: {total_downloads}"
    )

    for segment in SEGMENTS:
        print(
            f"  {segment}: {download_counts[segment]}"
        )

    print(
        f"Final alignment length: {final_length} nt"
    )

    print(
        f"Final MSA:\n{final_output}"
    )

    print(
        f"Pairwise identity matrix:\n{identity_matrix_path}"
    )

    print(
        f"Pairwise genetic-distance matrix:\n{distance_matrix_path}"
    )

    print(
        f"Within/between clade summary:\n{clade_summary_path}"
    )

    print(
        "\nOriginal *_NU.fasta files were not modified. "
        "Any sequences downloaded from NCBI were saved separately under "
        "00_downloaded_missing_sequences/."
    )


if __name__ == "__main__":
    main()
