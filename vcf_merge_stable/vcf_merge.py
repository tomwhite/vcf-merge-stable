from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from itertools import groupby
from typing import Iterator

import click
import cyvcf2

from .list_merge import merge_with
from .merge_none import can_merge, merge_record


@dataclass
class _Record:
    chrom: str
    pos: int
    id: str
    ref: str
    alts: list[str]
    qual: float | None
    filters: list[str]


def _from_variant(v: cyvcf2.Variant) -> _Record:
    return _Record(
        chrom=v.CHROM,
        pos=v.POS,
        id=v.ID or ".",
        ref=v.REF,
        alts=list(v.ALT) if v.ALT else ["."],
        qual=v.QUAL,
        filters=v.FILTER.split(";") if v.FILTER and v.FILTER != "PASS" else [],
    )


def _combine(a: _Record, b: _Record) -> _Record:
    _, alts = merge_record([(a.ref, a.alts), (b.ref, b.alts)])
    ids = list(dict.fromkeys(x for x in [a.id, b.id] if x != "."))
    quals = [q for q in [a.qual, b.qual] if q is not None]
    return _Record(
        chrom=a.chrom,
        pos=a.pos,
        id=";".join(ids) if ids else ".",
        ref=a.ref,
        alts=alts,
        qual=max(quals) if quals else None,
        filters=sorted(set(a.filters) | set(b.filters)),
    )


def _iter_positions(vcf: cyvcf2.VCF) -> Iterator[list[_Record]]:
    for _, group in groupby(vcf, key=lambda v: (v.CHROM, v.POS)):
        yield [_from_variant(v) for v in group]


def _merge_streams(
    s1: Iterator[list[_Record]],
    s2: Iterator[list[_Record]],
    contig_rank: dict[str, int],
) -> Iterator[tuple[list[_Record], list[_Record]]]:
    def key(recs: list[_Record]) -> tuple[int, int]:
        return contig_rank.get(recs[0].chrom, len(contig_rank)), recs[0].pos

    g1, g2 = next(s1, None), next(s2, None)
    while g1 is not None or g2 is not None:
        if g1 is None:
            yield [], g2
            g2 = next(s2, None)
        elif g2 is None:
            yield g1, []
            g1 = next(s1, None)
        elif key(g1) < key(g2):
            yield g1, []
            g1 = next(s1, None)
        elif key(g1) > key(g2):
            yield [], g2
            g2 = next(s2, None)
        else:
            yield g1, g2
            g1, g2 = next(s1, None), next(s2, None)


def _format_record(r: _Record) -> str:
    qual = f"{r.qual:.6g}" if r.qual is not None else "."
    filt = ";".join(r.filters) if r.filters else "PASS"
    return f"{r.chrom}\t{r.pos}\t{r.id}\t{r.ref}\t{','.join(r.alts)}\t{qual}\t{filt}\t."


def _parse_contigs(raw_header: str) -> list[tuple[str, str | None]]:
    result = []
    for line in raw_header.splitlines():
        m = re.match(r"##contig=<ID=([^,>]+)(?:,length=(\d+))?", line)
        if m:
            result.append((m.group(1), m.group(2)))
    return result


def merge_vcf_files(path1: str, path2: str, output: str | None = None) -> None:
    """Merge two VCF files, discarding INFO and FORMAT fields.

    Records at the same position are merged using -m none semantics (see merge_pairwise).
    Output contains only fixed fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO=.).
    When two records combine: IDs are joined with ";", QUAL takes the max, FILTER is the union.
    """
    vcf1 = cyvcf2.VCF(path1)
    vcf2 = cyvcf2.VCF(path2)

    vcf1_contig_ids: set[str] = set()
    all_contigs: list[tuple[str, str | None]] = []
    for cid, length in _parse_contigs(vcf1.raw_header):
        vcf1_contig_ids.add(cid)
        all_contigs.append((cid, length))

    extra_contigs: list[tuple[str, str | None]] = []
    vcf2_seen: set[str] = set()
    for cid, length in _parse_contigs(vcf2.raw_header):
        if cid not in vcf1_contig_ids and cid not in vcf2_seen:
            extra_contigs.append((cid, length))
            all_contigs.append((cid, length))
            vcf2_seen.add(cid)

    contig_rank = {cid: i for i, (cid, _) in enumerate(all_contigs)}

    with (open(output, "w") if output else sys.stdout) as out:
        for line in vcf1.raw_header.splitlines():
            if line.startswith("#CHROM"):
                for cid, length in extra_contigs:
                    out.write(f"##contig=<ID={cid}" + (f",length={length}" if length else "") + ">\n")
                fixed = "\t".join(line.split("\t")[:8])
                out.write(fixed + "\n")
                break
            out.write(line + "\n")

        for l1, l2 in _merge_streams(_iter_positions(vcf1), _iter_positions(vcf2), contig_rank):
            for r in merge_with(
                l1, l2,
                equiv=lambda a, b: can_merge(a.ref, a.alts, b.ref, b.alts),
                combine=_combine,
            ):
                out.write(_format_record(r) + "\n")


@click.command()
@click.argument("vcf1", type=click.Path(exists=True))
@click.argument("vcf2", type=click.Path(exists=True))
@click.option("-o", "--output", type=click.Path(), default=None, help="Output VCF file (default: stdout)")
def cli(vcf1: str, vcf2: str, output: str | None) -> None:
    """Merge two VCF files using -m none semantics with stable record ordering.

    Unlike bcftools merge -m none, output record order is stable: relative ordering
    from both input files is preserved (interleaved, not file-first). Records at the
    same position are merged when their alt sets overlap; all INFO and FORMAT fields
    are discarded.
    """
    merge_vcf_files(vcf1, vcf2, output)


if __name__ == "__main__":
    cli()
