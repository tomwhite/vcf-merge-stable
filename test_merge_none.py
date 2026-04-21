import os
import subprocess
import tempfile

import pytest

from merge_none import can_merge, group_records

BCFTOOLS = "bcftools"
BGZIP = "bgzip"
TABIX = "tabix"

VCF_TEMPLATE = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\t{ref}\t{alt}\t.\tPASS\t.
"""


def bcftools_merge_none(ref1: str, alt1: list[str], ref2: str, alt2: list[str]) -> bool:
    with tempfile.TemporaryDirectory() as d:
        for i, (ref, alt) in enumerate([(ref1, alt1), (ref2, alt2)], 1):
            vcf = os.path.join(d, f"s{i}.vcf")
            with open(vcf, "w") as f:
                f.write(VCF_TEMPLATE.format(ref=ref, alt=",".join(alt)))
            subprocess.run([BGZIP, vcf], check=True)
            subprocess.run([TABIX, "-p", "vcf", vcf + ".gz"], check=True)

        result = subprocess.run(
            [BCFTOOLS, "merge", "-m", "none",
             os.path.join(d, "s1.vcf.gz"),
             os.path.join(d, "s2.vcf.gz")],
            capture_output=True, text=True, check=True,
        )
        data_lines = [l for l in result.stdout.splitlines() if not l.startswith("#")]
        return len(data_lines) == 1


@pytest.mark.parametrize("ref1,alt1,ref2,alt2,expected", [
    # SNP cases
    ("A", ["T"], "A", ["C"], False),              # disjoint SNPs
    ("A", ["T", "C"], "A", ["C"], True),          # superset alts
    ("A", ["T"], "A", ["T"], True),               # identical records
    ("A", ["T"], "C", ["T"], False),              # different ref
    ("A", ["C"], "A", ["T", "C"], True),          # subset, either direction
    ("A", ["T", "G"], "A", ["C", "G"], True),     # shared allele, neither a subset
    ("A", ["T", "G"], "A", ["C"], False),         # disjoint multiallelic
    ("A", ["T", "C"], "A", ["T", "C"], True),     # identical multiallelics
    ("A", ["T", "C", "G"], "A", ["T", "C"], True),  # 3-allele superset
    # ref-only cases
    ("A", ["."], "A", ["."], True),               # both ref-only
    ("A", ["."], "A", ["T"], True),               # ref-only + SNP
    ("A", ["T"], "A", ["."], True),               # SNP + ref-only
    ("A", ["."], "A", ["T", "C"], True),          # ref-only + multiallelic
    # MNP cases
    ("AA", ["TT"], "AA", ["TT"], True),           # same MNP
    ("AA", ["TT"], "AA", ["CC"], False),          # different MNPs
    ("AA", ["TT", "CC"], "AA", ["CC"], True),     # multi MNP, shared
    ("AA", ["TT", "CC"], "AA", ["GG"], False),    # multi MNP, disjoint
    ("AA", ["TT"], "AA", ["T"], False),           # MNP vs deletion (different alt len)
    ("AA", ["TT"], "AA", ["TTA"], False),         # MNP vs insertion
    # indel cases
    ("A", ["AT"], "A", ["AC"], False),            # different insertions
    ("A", ["AT"], "A", ["AT"], True),             # same insertion
    ("AT", ["A"], "AT", ["A"], True),             # same deletion
    ("AT", ["A"], "AT", ["ATG"], False),          # deletion vs insertion, disjoint
    ("A", ["AT", "AC"], "A", ["AC"], True),       # multiallelic insertion, shared allele
    ("AT", ["A", "ATG"], "AT", ["ATG"], True),    # mixed indels, shared allele
    ("A", ["AT"], "AT", ["A"], False),            # different ref lengths
])
def test_can_merge(ref1, alt1, ref2, alt2, expected):
    assert can_merge(ref1, alt1, ref2, alt2) == expected


@pytest.mark.parametrize("ref1,alt1,ref2,alt2", [
    # SNP cases
    ("A", ["T"], "A", ["C"]),              # disjoint SNPs
    ("A", ["T", "C"], "A", ["C"]),         # superset alts
    ("A", ["T"], "A", ["T"]),              # identical records
    ("A", ["T"], "C", ["T"]),              # different ref
    ("A", ["C"], "A", ["T", "C"]),         # subset, either direction
    ("A", ["T", "G"], "A", ["C", "G"]),    # shared allele, neither a subset
    ("A", ["T", "G"], "A", ["C"]),         # disjoint multiallelic
    ("A", ["T", "C"], "A", ["T", "C"]),    # identical multiallelics
    ("A", ["T", "C", "G"], "A", ["T", "C"]),  # 3-allele superset
    # ref-only cases
    ("A", ["."], "A", ["."]),            # both ref-only
    ("A", ["."], "A", ["T"]),            # ref-only + SNP
    ("A", ["T"], "A", ["."]),            # SNP + ref-only
    ("A", ["."], "A", ["T", "C"]),       # ref-only + multiallelic
    # MNP cases
    ("AA", ["TT"], "AA", ["TT"]),          # same MNP
    ("AA", ["TT"], "AA", ["CC"]),          # different MNPs
    ("AA", ["TT", "CC"], "AA", ["CC"]),    # multi MNP, shared
    ("AA", ["TT", "CC"], "AA", ["GG"]),    # multi MNP, disjoint
    ("AA", ["TT"], "AA", ["T"]),           # MNP vs deletion (different alt len)
    ("AA", ["TT"], "AA", ["TTA"]),         # MNP vs insertion
    # indel cases
    ("A", ["AT"], "A", ["AC"]),           # different insertions
    ("A", ["AT"], "A", ["AT"]),           # same insertion
    ("AT", ["A"], "AT", ["A"]),           # same deletion
    ("AT", ["A"], "AT", ["ATG"]),         # deletion vs insertion, disjoint
    ("A", ["AT", "AC"], "A", ["AC"]),     # multiallelic insertion, shared allele
    ("AT", ["A", "ATG"], "AT", ["ATG"]),  # mixed indels, shared allele
    ("A", ["AT"], "AT", ["A"]),           # different ref lengths
])
def test_matches_bcftools(ref1, alt1, ref2, alt2):
    assert can_merge(ref1, alt1, ref2, alt2) == bcftools_merge_none(ref1, alt1, ref2, alt2)


# --- group_records tests ---

HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
)
ROW = "chr1\t100\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"


def bcftools_merge_none_n(
    file_records: list[list[tuple[str, list[str]]]]
) -> list[list[int]]:
    """Run bcftools merge -m none on N files (each a list of records) and return groups.

    Returns groups as lists of flat indices into the concatenation of all file_records.
    """
    with tempfile.TemporaryDirectory() as d:
        paths = []
        for i, recs in enumerate(file_records, 1):
            p = os.path.join(d, f"s{i}.vcf")
            with open(p, "w") as f:
                f.write(HDR.format(sample=f"S{i}"))
                for ref, alts in recs:
                    gt = "0/0" if alts == ["."] else "0/1"
                    f.write(ROW.format(ref=ref, alt=",".join(alts), gt=gt))
            subprocess.run([BGZIP, p], check=True)
            subprocess.run([TABIX, "-p", "vcf", p + ".gz"], check=True)
            paths.append(p + ".gz")

        result = subprocess.run(
            [BCFTOOLS, "merge", "-m", "none"] + paths,
            capture_output=True, text=True, check=True,
        )

    # Map each sample to its flat record indices (in file order)
    flat = [(ref, alts) for recs in file_records for ref, alts in recs]
    sample_indices: dict[str, list[int]] = {}
    offset = 0
    for i, recs in enumerate(file_records, 1):
        sample_indices[f"S{i}"] = list(range(offset, offset + len(recs)))
        offset += len(recs)

    header_line = next(l for l in result.stdout.splitlines() if l.startswith("#CHROM"))
    sample_cols = header_line.split("\t")[9:]

    # For each output record, find which input record each sample contributed.
    # Use the output record's ALT set to match: find the first unassigned input
    # record for each sample whose alts are a subset of the output alts.
    assigned: set[int] = set()
    groups: list[list[int]] = []
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        output_alts = set(fields[4].split(","))
        gts = dict(zip(sample_cols, fields[9:]))
        group = []
        for sample, gt in gts.items():
            if gt == "./.":
                continue
            for idx in sample_indices[sample]:
                if idx not in assigned:
                    _, alts = flat[idx]
                    if alts == ["."] or set(alts) <= output_alts:
                        group.append(idx)
                        assigned.add(idx)
                        break
        if group:
            groups.append(sorted(group))

    return groups


@pytest.mark.parametrize("records,expected_groups", [
    # two identical records
    ([("A", ["T"]), ("A", ["T"])],                        [[0, 1]]),
    # two disjoint records
    ([("A", ["T"]), ("A", ["C"])],                        [[0], [1]]),
    # three records: F1 anchor T pulls in F3, F2 alone
    ([("A", ["T"]), ("A", ["C"]), ("A", ["T", "C"])],     [[0, 2], [1]]),
    # chain: all share pairwise but anchor limits merging
    ([("A", ["T", "C"]), ("A", ["C", "G"]), ("A", ["G", "T"])],  [[0, 2], [1]]),
    # alt order in F1 changes anchor: C,T anchor=C pulls in F3 (G,C)
    ([("A", ["C", "T"]), ("A", ["G", "T"]), ("A", ["G", "C"])],  [[0, 2], [1]]),
    # all same: all merge
    ([("A", ["T"]), ("A", ["T"]), ("A", ["T"])],          [[0, 1, 2]]),
    # all disjoint: three separate groups
    ([("A", ["T"]), ("A", ["C"]), ("A", ["G"])],          [[0], [1], [2]]),
    # 2-file 3-record: F1 has two records, F2 has one matching the first
    ([("A", ["T"]), ("A", ["C"]), ("A", ["T"])],          [[0, 2], [1]]),
    # F2 multiallelic T,C: anchor T pulls in F1:[T], F1:[C] alone
    ([("A", ["T"]), ("A", ["C"]), ("A", ["T", "C"])],     [[0, 2], [1]]),
    # ref-only merges into current group
    ([("A", ["T"]), ("A", ["."]), ("A", ["C"])],          [[0, 1], [2]]),
])
def test_group_records(records, expected_groups):
    assert group_records(records) == expected_groups


@pytest.mark.parametrize("file_records", [
    # 3 records across 2 files
    [[("A", ["T"]), ("A", ["C"])], [("A", ["T"])]],
    [[("A", ["T"]), ("A", ["C"])], [("A", ["C"])]],
    [[("A", ["T"]), ("A", ["C"])], [("A", ["G"])]],
    [[("A", ["T"]), ("A", ["C"])], [("A", ["T", "C"])]],
    [[("A", ["T", "C"])],          [("A", ["T"]), ("A", ["C"])]],
    # 3 records across 3 files
    [[("A", ["T"])], [("A", ["C"])], [("A", ["T", "C"])]],
    [[("A", ["T", "C"])], [("A", ["C", "G"])], [("A", ["G", "T"])]],
    [[("A", ["C", "T"])], [("A", ["G", "T"])], [("A", ["G", "C"])]],
])
def test_group_records_matches_bcftools(file_records):
    flat = [(ref, alts) for recs in file_records for ref, alts in recs]
    assert {frozenset(g) for g in group_records(flat)} == {frozenset(g) for g in bcftools_merge_none_n(file_records)}
