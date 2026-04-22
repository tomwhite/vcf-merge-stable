# vcf-merge-stable

Merge two sample-less VCF files using `bcftools merge -m none` semantics, with one
important difference: **output record order is stable**.

## The problem

`bcftools merge -m none` does not preserve the relative ordering of input records.
Consider two files with three distinct ALT alleles at the same position:

**f1.vcf** — A/T then A/C  
**f2.vcf** — A/G then A/C

```
# f1.vcf               # f2.vcf
chr1  100  A  T        chr1  100  A  G
chr1  100  A  C        chr1  100  A  C
```

Running `bcftools merge -m none example/f1.vcf.gz example/f2.vcf.gz`:

```
chr1  100  A  C   ← C appears first, inconsistent with both input files
chr1  100  A  T
chr1  100  A  G
```

The output order is consistent with neither file. A/C appears before A/T even though
A/T is first in f1, and A/C appears before A/G even though A/G is first in f2.

## The fix

`vcf-merge-stable` uses a topological sort to interleave records from both files,
preserving relative ordering from each:

```
$ vcf-merge-stable example/f1.vcf.gz example/f2.vcf.gz

chr1  100  A  T   ← first in f1
chr1  100  A  G   ← first in f2
chr1  100  A  C   ← second in both (merged into one record)
```

The output respects the ordering constraints from both inputs: T before C (from f1),
G before C (from f2). A/C appears exactly once because the identical records from
each file are merged.

## Usage

```
vcf-merge-stable [OPTIONS] VCF1 VCF2

Options:
  -o, --output PATH  Output VCF file (default: stdout)
```

Input files must be bgzipped and tabix-indexed. Output contains only fixed VCF fields
(CHROM, POS, ID, REF, ALT, QUAL, FILTER); INFO and FORMAT are discarded.

When two records are merged: IDs are joined with `;`, QUAL takes the max, FILTER is
the union.

## Mergeability

Two records at the same position are merged when:

1. They have the same REF allele
2. Their ALT sets have a non-empty intersection (share at least one alt allele)
3. A ref-only site (`ALT=.`) merges with anything

This matches `bcftools merge -m none` pairwise mergeability exactly.
