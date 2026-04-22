from .list_merge import merge_with


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


def merge_pairwise(
    l1: list[tuple[str, list[str]]],
    l2: list[tuple[str, list[str]]],
) -> list[tuple[str, list[str]]]:
    """Merge records from two files at the same position.

    Uses merge_with to preserve relative ordering from both files, matching
    records via can_merge and combining matched pairs via merge_record.
    For N files, fold: functools.reduce(merge_pairwise, file_record_lists).
    """
    return merge_with(
        l1, l2,
        equiv=lambda a, b: can_merge(a[0], a[1], b[0], b[1]),
        combine=lambda a, b: merge_record([a, b]),
    )
