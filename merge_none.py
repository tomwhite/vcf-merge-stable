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
