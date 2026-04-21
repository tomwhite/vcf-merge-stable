import pytest

from merge_none import merge_alleles, merge_record, merge_records


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


@pytest.mark.parametrize("records,expected", [
    # all identical: one output record
    ([("A", ["T"]), ("A", ["T"])],                          [("A", ["T"])]),
    # disjoint: two separate records preserved
    ([("A", ["T"]), ("A", ["C"])],                          [("A", ["T"]), ("A", ["C"])]),
    # anchor T pulls in T,C; C alone
    ([("A", ["T"]), ("A", ["C"]), ("A", ["T", "C"])],       [("A", ["T", "C"]), ("A", ["C"])]),
    # chain: anchor T pulls in G,T; C,G alone
    ([("A", ["T", "C"]), ("A", ["C", "G"]), ("A", ["G", "T"])], [("A", ["T", "C", "G"]), ("A", ["C", "G"])]),
    # ref-only merges into first group
    ([("A", ["T"]), ("A", ["."]), ("A", ["C"])],            [("A", ["T"]), ("A", ["C"])]),
    # all ref-only
    ([("A", ["."]), ("A", ["."])],                          [("A", ["."])]),
    # single record passthrough
    ([("A", ["T"])],                                        [("A", ["T"])]),
])
def test_merge_records(records, expected):
    assert merge_records(records) == expected
