import os
import subprocess
import tempfile

import pytest

from merge_none import merge_alleles, merge_record, merge_pairwise

BCFTOOLS = "bcftools"
BGZIP = "bgzip"
TABIX = "tabix"

VCF_HDR = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)
VCF_ROW = "chr1\t100\t.\t{ref}\t{alt}\t.\tPASS\t.\n"


def bcftools_merge_pairwise(
    l1: list[tuple[str, list[str]]],
    l2: list[tuple[str, list[str]]],
) -> set[tuple[str, frozenset[str]]]:
    """Run bcftools merge -m none on two sample-less files, return output records as a set."""
    with tempfile.TemporaryDirectory() as d:
        paths = []
        for i, recs in enumerate([l1, l2], 1):
            p = os.path.join(d, f"s{i}.vcf")
            with open(p, "w") as f:
                f.write(VCF_HDR)
                for ref, alts in recs:
                    f.write(VCF_ROW.format(ref=ref, alt=",".join(alts)))
            subprocess.run([BGZIP, p], check=True)
            subprocess.run([TABIX, "-p", "vcf", p + ".gz"], check=True)
            paths.append(p + ".gz")

        result = subprocess.run(
            [BCFTOOLS, "merge", "-m", "none"] + paths,
            capture_output=True, text=True, check=True,
        )

    return {
        (fields[3], frozenset(fields[4].split(",")))
        for line in result.stdout.splitlines()
        if not line.startswith("#")
        for fields in [line.split("\t")]
    }


@pytest.mark.parametrize("alt_lists,expected", [
    # single record
    ([["T"]],               ["T"]),           # single biallelic
    ([["T", "C"]],          ["T", "C"]),      # single multiallelic
    ([["."]], ["."]  ),                        # single ref-only
    # identical records
    ([["T"], ["T"]],        ["T"]),            # same allele deduplicated
    ([["T", "C"], ["T", "C"]], ["T", "C"]),   # same multiallelic deduplicated
    # disjoint (shouldn't normally be grouped, but merge_alleles is not responsible for that)
    ([["T"], ["C"]],        ["T", "C"]),       # disjoint: order from first record first
    # overlapping
    ([["T", "G"], ["C", "G"]], ["T", "G", "C"]),  # shared G, new C appended
    ([["T", "C"], ["C", "G"]], ["T", "C", "G"]),  # shared C, new G appended
    ([["T"], ["T", "C"]],   ["T", "C"]),       # subset: C appended
    ([["T", "C"], ["T"]],   ["T", "C"]),       # superset: nothing new
    # ref-only mixed in
    ([["T"], ["."]],        ["T"]),            # ref-only contributes nothing
    ([["."], ["T"]],        ["T"]),            # ref-only first, then real allele
    ([["T"], ["."], ["C"]], ["T", "C"]),       # ref-only in the middle
    ([["."], ["."]],        ["."]),            # all ref-only
    # three records
    ([["T", "C"], ["C", "G"], ["G", "T"]], ["T", "C", "G"]),  # chain
    ([["T"], ["C"], ["G"]], ["T", "C", "G"]),                  # all disjoint
])
def test_merge_alleles(alt_lists, expected):
    assert merge_alleles(alt_lists) == expected


@pytest.mark.parametrize("records,expected", [
    # single record passthrough
    ([("A", ["T"])],                        ("A", ["T"])),
    ([("A", ["T", "C"])],                   ("A", ["T", "C"])),
    ([("A", ["."])],                        ("A", ["."])),
    # identical records
    ([("A", ["T"]), ("A", ["T"])],          ("A", ["T"])),
    # overlapping
    ([("A", ["T", "G"]), ("A", ["C", "G"])], ("A", ["T", "G", "C"])),
    ([("A", ["T"]), ("A", ["T", "C"])],     ("A", ["T", "C"])),
    # ref-only mixed in
    ([("A", ["T"]), ("A", ["."])],          ("A", ["T"])),
    ([[("A", ["."])][0], ("A", ["T"])],     ("A", ["T"])),
    ([("A", ["."]), ("A", ["."])],          ("A", ["."])),
    # indels
    ([("AT", ["A"]), ("AT", ["A"])],        ("AT", ["A"])),
    ([("A", ["AT"]), ("A", ["AT", "AC"])],  ("A", ["AT", "AC"])),
])
def test_merge_record(records, expected):
    assert merge_record(records) == expected


@pytest.mark.parametrize("l1,l2,expected", [
    # no matches: all records separate, ordering from both files preserved
    ([("A", ["T"])],        [("A", ["C"])],              [("A", ["T"]), ("A", ["C"])]),
    # identical records merge
    ([("A", ["T"])],        [("A", ["T"])],              [("A", ["T"])]),
    # shared allele: records merge
    ([("A", ["T", "G"])],   [("A", ["C", "G"])],         [("A", ["T", "G", "C"])]),
    # ref-only merges with anything
    ([("A", ["T"])],        [("A", ["."])],              [("A", ["T"])]),
    # two records each, one match: ordering T, G, C (not T, C, G)
    ([("A", ["T"]), ("A", ["C"])], [("A", ["G"]), ("A", ["C"])],
     [("A", ["T"]), ("A", ["G"]), ("A", ["C"])]),
    # l1 multiallelic matches first l2 record; second l2 record separate
    ([("A", ["T", "C"])],   [("A", ["T"]), ("A", ["C"])],
     [("A", ["T", "C"]), ("A", ["C"])]),
    # l1 has two records, l2 multiallelic matches first
    ([("A", ["T"]), ("A", ["C"])], [("A", ["T", "C"])],
     [("A", ["T", "C"]), ("A", ["C"])]),
    # different ref: no merge
    ([("A", ["T"])],        [("C", ["T"])],              [("A", ["T"]), ("C", ["T"])]),
])
def test_merge_pairwise(l1, l2, expected):
    assert merge_pairwise(l1, l2) == expected


@pytest.mark.parametrize("l1,l2", [
    ([("A", ["T"])],        [("A", ["C"])]),
    ([("A", ["T"])],        [("A", ["T"])]),
    ([("A", ["T", "G"])],   [("A", ["C", "G"])]),
    ([("A", ["T"])],        [("A", ["."])]),
    ([("A", ["T"]), ("A", ["C"])], [("A", ["G"]), ("A", ["C"])]),
    ([("A", ["T", "C"])],   [("A", ["T"]), ("A", ["C"])]),
    ([("A", ["T"]), ("A", ["C"])], [("A", ["T", "C"])]),
])
def test_merge_pairwise_matches_bcftools(l1, l2):
    expected = bcftools_merge_pairwise(l1, l2)
    actual = {(ref, frozenset(alts)) for ref, alts in merge_pairwise(l1, l2)}
    assert actual == expected
