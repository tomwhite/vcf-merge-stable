def merge_none(ref1: str, alt1: list[str], ref2: str, alt2: list[str]) -> bool:
    """Return True if two VCF records at the same position can be merged under bcftools merge -m none.

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


def merge_none_n(records: list[tuple[str, list[str]]]) -> list[list[int]]:
    """Group VCF records at the same position into merge groups under bcftools merge -m none.

    Takes a flat ordered list of (ref, alts) records — as bcftools sees them, ordered
    by input file then by position within that file — and returns a list of groups,
    where each group is a list of record indices that bcftools would merge into one
    output record.

    The algorithm mirrors bcftools' greedy sequential approach:
    1. Take the first alt allele of the first unprocessed record as the anchor.
    2. Collect all subsequent unprocessed records whose alts contain the anchor
       (or are ref-only).
    3. Emit that set as one group, then repeat from step 1.

    A ref-only record (alts == ["."]) is absorbed into the first group whose anchor
    it could match — i.e. it joins the current group being built.
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
