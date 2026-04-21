import pytest

from merge_none import merge_alleles


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
