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

Case 4 is the clearest illustration: `F2:[T,C]` anchors on T (its first alt), pulls in
`F1:[T]`, and leaves `F1:[C]` alone — even though C is also in the multiallelic.

## Extending the Python implementation

Extending to N records requires this anchor-based logic, treating all records at the
same position as a flat list regardless of which file they come from. A suitable
signature would be:

```python
def merge_none_n(records: list[tuple[str, list[str]]]) -> list[list[int]]:
    ...  # returns groups of record indices
```
