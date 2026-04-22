import os
import subprocess
import tempfile

import pytest

from vcf_merge_stable.merge_none import can_merge

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
