Rotavirus P[8]DS-1 Signature and Allele Analysis

This repository contains the sequence files, analysis scripts, and tabular outputs used to evaluate allele structure, clade relationships, genetic signatures, allele enrichment, and nucleotide genetic distances in emergent P[8]DS-1 human rotaviruses.

The analyses focus on VP4, VP6, VP1, VP2, VP3, NSP1, NSP2, NSP3, and NSP4. VP7 and NSP5/6 were not included in the allele-signature analyses described here.

Important file-naming note

Two different sets of nucleotide FASTA files are used for different analyses. They should not be treated as interchangeable.

Files such as VP1_NU.fasta in the clade analysis dataset are used for the concatenated clade-level analysis.

Files such as VP1_NU (1).fasta were used for genetic-signature detection.

For GitHub, it is strongly recommended to rename these files so their purpose is immediately clear, for example:

VP1_NU_clade.fasta
VP1_NU_signature.fasta

The same naming convention can be applied to the other segments.

Repository organization

A suggested organization is:

Rotavirus-P8DS1-Signature-Analysis/
│
├── README.md
├── requirements.txt
│
├── clade_analysis/
│   ├── VP4_NU_clade.fasta
│   ├── VP6_NU_clade.fasta
│   ├── VP1_NU_clade.fasta
│   ├── VP2_NU_clade.fasta
│   ├── VP3_NU_clade.fasta
│   ├── NSP1_NU_clade.fasta
│   ├── NSP2_NU_clade.fasta
│   ├── NSP3_NU_clade.fasta
│   ├── NSP4_NU_clade.fasta
│   ├── Clades_accessions.csv
│   ├── P1_build_concatenated_MSA_with_clade_distances.py
│   └── outputs/
│
├── genetic_signature_analysis/
│   ├── VP1_AA.fasta
│   ├── VP1_NU_signature.fasta
│   ├── ...
│   ├── VP4_AA.fasta
│   ├── VP4_NU_signature.fasta
│   ├── ...
│   ├── NSP4_AA.fasta
│   ├── NSP4_NU_signature.fasta
│   ├── P1_signature_discovery_callable_only.py
│   └── outputs/
│
├── allele_analysis/
│   ├── strain_group_allele_assignments.csv
│   └── related allele-assignment tables
│
├── fisher_test/
│   ├── VP4_fisher_test.csv
│   ├── VP1_NSP4_fisher_test.csv
│   ├── Fisher_test_Table_S4_manuscript_ready.py
│   └── Table_S4_Fisher_allele_enrichment.csv
│
└── genetic_distance_analysis/
    ├── VP1.csv
    ├── VP2.csv
    ├── VP3.csv
    ├── VP4.csv
    ├── VP6.csv
    ├── NSP1.csv
    ├── NSP2.csv
    ├── NSP3.csv
    ├── NSP4.csv
    ├── P1_allele_genetic_distance_analysis.py
    └── outputs/

1. Clade analysis and concatenated genome alignment

Purpose

This analysis was used to quantitatively evaluate the genetic separation of P[8]DS-1 Clades I-V.

The nucleotide sequences for the following nine segments were analyzed:

VP4
VP6
VP1
VP2
VP3
NSP1
NSP2
NSP3
NSP4

Each segment was aligned separately using MAFFT and the aligned segments were concatenated strain-by-strain in the fixed order:

VP4–VP6–VP1–VP2–VP3–NSP1–NSP2–NSP3–NSP4

The concatenated alignment was then used to calculate pairwise nucleotide identity and genetic distance among strains assigned to Clades I-V.

Genetic distance was defined as:

Genetic distance = 100 - nucleotide identity (%)

Mean within-clade and between-clade nucleotide identities and genetic distances were then compared.

Main inputs

The nucleotide FASTA files in this section are the files used specifically for the clade/concatenation analysis.

Example:

VP1_NU_clade.fasta

This file is distinct from the VP1 nucleotide FASTA used for genetic-signature detection.

Main outputs

Depending on the version of the script used, outputs include:

P1_reassortant_concatenated_MSA.fasta
P1_concatenated_MSA_partitions.csv
P1_concatenated_MSA_sequence_sources.csv
P1_clades_pairwise_identity_matrix.csv
P1_clades_pairwise_genetic_distance_matrix.csv
P1_clade_within_between_identity_distance_summary.csv

The final within/between summary was used to quantitatively support the Clades I-V assignments.

2. Genetic-signature detection

Purpose

This analysis identifies nucleotide and amino-acid changes associated with type-B alleles and the P[8]DS-1 reassortant background.

The analysis includes:

VP1
VP2
VP3
VP4
VP6
NSP1
NSP2
NSP3
NSP4

Both nucleotide (NU) and amino-acid (AA) alignments are used.

The files originally named with duplicate suffixes such as:

VP1_NU (1).fasta

belong to this genetic-signature analysis. To avoid confusion in the repository, these should preferably be renamed:

VP1_NU_signature.fasta

FASTA header groups

For non-VP4 segments:

RP_   = P[8]DS-1 reassortant carrying the type-B allele
WTP_  = control/wild-type strain carrying the type-B allele
R_    = type-A allele
O_    = type-C allele
G_    = type-D allele

For VP4:

LGR_  = P[8]DS-1 reassortant / light-green P[8]-B group
LGWT_ = control strain carrying the light-green P[8]-B-like allele
DG_   = deep-green P[8]-A group

Signature criteria

Criterion 1

For non-VP4 genes, the candidate state must be:

100% in RP_
100% in WTP_
0% in R_, O_, and G_ where those groups are present

For VP4:

100% in LGR_
100% in LGWT_
0% in DG_

Criterion 2

For non-VP4 genes:

100% in RP_
>0% but <100% in WTP_
0% in R_, O_, and G_ where those groups are present

For VP4:

100% in LGR_
>0% but <100% in LGWT_
0% in DG_

Gaps and ambiguous characters are treated as non-callable. Percentages are calculated using callable sequences as the denominator.

Main script

P1_signature_discovery_callable_only.py

Main outputs

P1_AA_signature_changes.csv
P1_NU_signature_changes.csv

3. Allele assignment analysis

Purpose

Intra-genotypic allele groups were defined from the phylogenetic clustering of each gene and represented as types A-D.

The allele designations used throughout the downstream analyses are:

A
B
C
D

For the original color-coded allele tables:

red    = A
pink   = B
orange = C
grey   = D

For VP4:

deep green  = P[8]-A
light green = P[8]-B

The allele-assignment tables link each strain to its viral population and its allele assignment for each analyzed segment.

Example columns:

Strain
Group
VP4
VP6
VP1
VP2
VP3
NSP1
NSP2
NSP3
NSP4

These allele assignments were subsequently used for the Fisher exact tests and allele-frequency analyses.

4. Per-gene Fisher exact tests

Purpose

Per-gene Fisher exact tests were used to determine whether the observed allele distributions differed significantly between P[8]DS-1 reassortants and the corresponding control viral populations.

Only strains with an assigned allele for the corresponding segment are included. Blank or unassigned entries are excluded from the denominator.

VP4

VP4 is analyzed separately because both groups contain genotype P[8] VP4 genes.

The comparison is:

P[8]-B versus P[8]-A

between:

P[8]DS-1
and
P[8]Wa

The 2 x 2 Fisher table therefore evaluates the distribution of P[8]-A and P[8]-B between these two viral populations.

VP6, VP1-VP3, and NSP1-NSP4

For these genotype-2 backbone segments, type-B alleles are compared with the non-B allele classes actually observed for that segment.

Examples:

VP6   : B vs A
VP1   : B vs A+C+D
VP2   : B vs A
VP3   : B vs A+C
NSP1  : B vs A
NSP2  : B vs A+C
NSP3  : B vs A
NSP4  : B vs A

The viral populations compared are:

P[8]DS-1
versus
P[4]DS-1

Two-sided Fisher's exact tests are used.

P values are corrected across the nine gene-level comparisons using the Benjamini-Hochberg false-discovery-rate procedure.

Main inputs

VP4_fisher_test.csv
VP1_NSP4_fisher_test.csv

Main script

Fisher_test_Table_S4_manuscript_ready.py

Main output

Table_S4_Fisher_allele_enrichment.csv

The output reports:

Segment
Allele comparison
P[8]DS-1 n/N (%)
Control group
Control n/N (%)
Odds ratio
Raw Fisher P value
Benjamini-Hochberg adjusted P value
Significance before correction
Significance after correction

5. Allele-level nucleotide identity and genetic-distance analysis

Purpose

This analysis provides quantitative support for the phylogenetically defined A-D allele groups.

Pairwise nucleotide identity is calculated within and between allele groups for each segment.

Genetic distance is defined as:

Genetic distance = 100 - nucleotide identity (%)

The analysis evaluates:

within-B identity/distance
within-other-allele identity/distance
B-versus-other-allele identity/distance

where the comparison allele can be A, C, or D depending on the segment.

The analysis also reports:

Identity ratio = mean within-B identity / mean B-versus-other identity

Distance ratio = mean within-B distance / mean B-versus-other distance

Allele separation is supported when within-allele identity is higher than between-allele identity and within-allele genetic distance is lower than the corresponding between-allele distance.

Main inputs

Segment-level CSV files:

VP1.csv
VP2.csv
VP3.csv
VP4.csv
VP6.csv
NSP1.csv
NSP2.csv
NSP3.csv
NSP4.csv

Main script

P1_allele_genetic_distance_analysis.py

Main output

P1_allele_genetic_distance_summary.csv

Additional pairwise tables or plots may also be produced depending on the version of the analysis script.

Software requirements

The analyses were performed using Python 3 and standard scientific Python packages. Depending on the script, required packages may include:

numpy
scipy
matplotlib

MAFFT is required for the concatenated clade-analysis workflow.

Package versions used for reproducible execution should be recorded in:

requirements.txt

Reproducibility notes

The clade and signature FASTA files are separate datasets and should be clearly distinguished by filename or directory.

Sequence headers should not be altered unless the corresponding grouping logic in the scripts is also updated.

For signature detection, allele-group prefixes in the FASTA headers are required for sequence classification.

For Fisher analyses, blank/unassigned allele cells are excluded rather than treated as an allele class.

Benjamini-Hochberg correction is applied across the nine per-gene Fisher tests.

Some scripts developed on the Wake Forest University DEAC cluster contain local file paths. Users running the analysis elsewhere may need to update those paths.

Citation

If you use this repository, please cite the associated manuscript:

Otoo, L. L. and Esstman, S. M.
Genetic Signatures of Emergent DS-1-like Reassortant Human Rotaviruses Exhibiting Wa-like P[8] VP4 Genes.

Citation details will be updated following publication.
