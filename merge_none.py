from list_merge import merge_with


def can_merge(ref1: str, alt1: list[str], ref2: str, alt2: list[str]) -> bool:
    """Return True if two VCF records at the same position can be merged under bcftools -m none.

    The rule is non-empty intersection of alt allele sets: records can merge when they
    share at least one alt allele string. This is more permissive than a subset check —
    e.g. A/[T,G] and A/[C,G] merge because they share "G", producing a new triallelic
    A/[T,G,C].

    Records with completely disjoint alt sets (e.g. ["T"] and ["C"]) cannot merge;
    bcftools emits them as two separate records with missing genotypes for the other sample.

    A ref-only site (ALT=".") always merges with any other record at the same position,
    since it introduces no new alleles.

    Verified against bcftools 1.23.1 including indel cases: pure deletions, pure
    insertions, and mixed del+ins or ins+SNP multiallelics all follow the same
    non-empty intersection rule. (bcftools 1.20 had different behaviour for mixed-type
    multiallelics, but this was a bug that was fixed in later versions.)
    """
    if ref1 != ref2:
        return False
    if alt1 == ["."] or alt2 == ["."]:
        return True
    return bool(set(alt1) & set(alt2))


def merge_alleles(alt_lists: list[list[str]]) -> list[str]:
    """Return the merged ALT list for a group of records, in first-occurrence order.

    Iterates alt_lists in order, appending any allele not yet seen. Ref-only
    records (["."]) contribute no alleles. Returns ["."] if all records are ref-only.
    """
    seen: list[str] = []
    seen_set: set[str] = set()
    for alts in alt_lists:
        if alts == ["."]:
            continue
        for a in alts:
            if a not in seen_set:
                seen.append(a)
                seen_set.add(a)
    return seen if seen else ["."]


def merge_record(records: list[tuple[str, list[str]]]) -> tuple[str, list[str]]:
    """Merge a pre-grouped list of (ref, alts) records into a single (ref, alts) record.

    All records must share the same REF. ALT alleles are merged via merge_alleles.
    """
    ref = records[0][0]
    return ref, merge_alleles([alts for _, alts in records])


def group_records(records: list[tuple[str, list[str]]]) -> list[list[int]]:
    """Group VCF records at the same position into merge groups under bcftools -m none.

    Takes a flat ordered list of (ref, alts) records, ordered by a k-way merge of the
    input files (preserving relative order within each file), and returns a list of groups,
    where each group is a list of record indices that bcftools would merge into one
    output record.

    The algorithm mirrors bcftools' greedy sequential approach:
    1. Take the first alt allele of the first unprocessed record as the anchor.
    2. Collect all subsequent unprocessed records whose alts contain the anchor
       (or are ref-only).
    3. Emit that set as one group, then repeat from step 1.

    A ref-only record (alts == ["."]) is absorbed into the first group whose anchor
    it could match — i.e. it joins the current group being built.

    Output order is guaranteed to follow input order: groups appear in the order their
    anchors are encountered, and indices within each group are in ascending order.
    This differs from bcftools, which outputs groups in k-way merge order across files
    and does not preserve input record ordering.
    """
    remaining = list(range(len(records)))
    groups: list[list[int]] = []

    while remaining:
        anchor_idx = remaining[0]
        _, anchor_alts = records[anchor_idx]
        if anchor_alts == ["."]:
            anchor = None  # ref-only anchor: accepts everything
        else:
            anchor = anchor_alts[0]

        group = [anchor_idx]
        leftover = []
        for idx in remaining[1:]:
            ref, alts = records[idx]
            if anchor is None or alts == ["."] or anchor in alts:
                group.append(idx)
            else:
                leftover.append(idx)

        groups.append(group)
        remaining = leftover

    return groups


def merge_two_files(
    l1: list[tuple[str, list[str]]],
    l2: list[tuple[str, list[str]]],
) -> list[tuple[str, list[str]]]:
    """Merge records from two files at the same position.

    Uses merge_with to preserve relative ordering from both files, matching
    records via can_merge and combining matched pairs via merge_record.
    """
    return merge_with(
        l1, l2,
        equiv=lambda a, b: can_merge(a[0], a[1], b[0], b[1]),
        combine=lambda a, b: merge_record([a, b]),
    )


def merge_records(records: list[tuple[str, list[str]]]) -> list[tuple[str, list[str]]]:
    """Group and merge all records at the same position, returning one record per group.

    Implements bcftools -m none semantics: groups via the greedy anchor algorithm
    (see group_records), then merges each group's alleles (see merge_record).
    Output order follows input order of each group's anchor record.
    """
    return [merge_record([records[i] for i in group]) for group in group_records(records)]
