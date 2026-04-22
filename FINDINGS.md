# bcftools merge -m none: findings

## Two-record rule

Two records at the same position can be merged if:

1. They have the same REF allele
2. Their ALT sets have a non-empty intersection (share at least one alt allele string)
3. A ref-only site (`ALT=.`) merges with anything

The intersection rule is more permissive than a subset check: `A/[T,G]` and `A/[C,G]`
merge because they share `G`, producing a new triallelic `A/[T,G,C]`.

Verified against bcftools 1.23.1. (bcftools 1.20 had different behaviour for mixed
del+ins multiallelics, but this was a bug fixed in later versions.)

## N-record rule (3+ records at the same position)

bcftools does **not** use connected components of pairwise mergeability. Instead it
uses a greedy sequential algorithm anchored on each record's first alt allele:

1. Take the **first alt allele** of the first unprocessed record as the anchor
2. Scan subsequent unprocessed records; merge all that **contain** the anchor allele
3. Repeat from step 1 with the remaining unprocessed records

This is order-dependent and can give different results to connected components:

| Input | Pairwise CC predicts | bcftools produces |
|---|---|---|
| `T` + `C` + `T,C` | all 3 merge | F1+F3 (anchor T∈T,C), F2 alone |
| `T,C` + `C,G` + `G,T` | all 3 merge (chain) | F1+F3 (anchor T∈G,T), F2 alone |
| `C,T` + `G,T` + `G,C` | all 3 merge (chain) | F1+F3 (anchor C∈G,C), F2 alone |

Note: swapping the alt order in the first record (e.g. `T,C` → `C,T`) changes which
subsequent records get merged, because only the **first** alt is used as the anchor.

**Biallelic special case:** when all records are biallelic (single ALT each), the anchor
algorithm and connected components always agree. With one alt per record, "shares an alt"
means "has the same alt", so mergeability is just equality — which is transitive by
definition. Ordering effects and anchor/CC divergence only arise when multiallelics are
present.

## 2-file, 3-record case (one file contributes 2 records)

The anchor algorithm treats all records at the same position as a single flat list,
ordered by file then by position within the file. File boundaries don't matter.

| File 1 records | File 2 record | Output |
|---|---|---|
| `A/T`, `A/C` | `A/T` | `A/T` (S1+S2), `A/C` (S1 only) |
| `A/T`, `A/C` | `A/C` | `A/C` (S1+S2), `A/T` (S1 only) |
| `A/T`, `A/C` | `A/G` | 3 separate records |
| `A/T`, `A/C` | `A/T,C` | `A/T,C` (S1-T + S2), `A/C` (S1 only) |
| `A/T,C` | `A/T`, `A/C` | `A/T,C` (S1 + S2-T), `A/C` (S2 only) |

Case 4 is the clearest illustration: `F2:[A/T,C]` anchors on T (its first alt), pulls in
`F1:[A/T]`, and leaves `F1:[A/C]` alone — even though C is also in the multiallelic.

**Output ordering:** bcftools does not preserve input record order in its output. For
example, `F1:[A/T], F1:[A/C], F2:[A/C]` — the anchor algorithm groups `{A/T}` first and
`{A/C, A/C}` second, but bcftools emits `A/C` before `A/T` (k-way merge order across
files). The Python `group_records` deliberately preserves input order: group order and
within-group order both follow the original record sequence.

**Input ordering:** the flat list passed to `group_records` should be produced by a
k-way merge of the input files that preserves relative ordering within each file.
Simple file-first concatenation gives the wrong output order. For example,
`F1:[A/T, A/C]` and `F2:[A/G, A/C]` concatenated as `[A/T, A/C, A/G, A/C]` yields
`A/T, A/C, A/G` — inconsistent with file-2 record order. The correct merged order
`[A/T, A/G, A/C, A/C]` yields `A/T, A/G, A/C`. When records from different files
are at the same position they conflict (no natural ordering between them), so the
k-way merge needs a defined tiebreaking rule.

`merge_pairwise` avoids this problem entirely for the 2-file case: it takes two
lists directly and uses a topological sort (via `merge_with`) to derive a consistent
output order from both inputs, without requiring pre-interleaving.

## Python implementation

The following functions are implemented in `merge_none.py`:

- `can_merge(ref1, alt1, ref2, alt2) -> bool` — pairwise mergeability check
- `merge_alleles(alt_lists) -> list[str]` — union of ALT alleles in first-occurrence order
- `merge_record(records) -> (ref, alts)` — merge a pre-grouped list of records
- `merge_pairwise(l1, l2) -> list[(ref, alts)]` — merge records from two files,
  preserving relative ordering from both via topological sort; fold over files for N>2
