import os
import subprocess
import tempfile

import pytest

from vcf_merge_stable.vcf_merge import merge_vcf_files

BGZIP = "bgzip"
TABIX = "tabix"

HDR = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
ROW = "chr1\t100\t{id}\t{ref}\t{alt}\t{qual}\t{filt}\t.\n"


def make_vcf(d: str, name: str, rows: list[dict]) -> str:
    p = os.path.join(d, f"{name}.vcf")
    with open(p, "w") as f:
        f.write(HDR)
        for row in rows:
            f.write(ROW.format(
                id=row.get("id", "."),
                ref=row["ref"],
                alt=",".join(row["alts"]),
                qual=row.get("qual", "."),
                filt=row.get("filt", "PASS"),
            ))
    subprocess.run([BGZIP, p], check=True)
    subprocess.run([TABIX, "-p", "vcf", p + ".gz"], check=True)
    return p + ".gz"


def run_merge(rows1: list[dict], rows2: list[dict]) -> list[dict]:
    with tempfile.TemporaryDirectory() as d:
        p1 = make_vcf(d, "f1", rows1)
        p2 = make_vcf(d, "f2", rows2)
        out = os.path.join(d, "out.vcf")
        merge_vcf_files(p1, p2, out)
        records = []
        with open(out) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                records.append({
                    "ref": fields[3],
                    "alts": fields[4].split(","),
                    "qual": fields[5],
                    "filt": fields[6],
                    "id": fields[2],
                })
        return records


def test_disjoint_records():
    result = run_merge(
        [{"ref": "A", "alts": ["T"]}],
        [{"ref": "A", "alts": ["C"]}],
    )
    assert len(result) == 2
    assert result[0]["alts"] == ["T"]
    assert result[1]["alts"] == ["C"]


def test_identical_records_merge():
    result = run_merge(
        [{"ref": "A", "alts": ["T"]}],
        [{"ref": "A", "alts": ["T"]}],
    )
    assert result == [{"ref": "A", "alts": ["T"], "qual": ".", "filt": "PASS", "id": "."}]


def test_overlapping_records_merge():
    result = run_merge(
        [{"ref": "A", "alts": ["T", "G"]}],
        [{"ref": "A", "alts": ["C", "G"]}],
    )
    assert len(result) == 1
    assert result[0]["ref"] == "A"
    assert result[0]["alts"] == ["T", "G", "C"]


def test_qual_max():
    result = run_merge(
        [{"ref": "A", "alts": ["T"], "qual": "30"}],
        [{"ref": "A", "alts": ["T"], "qual": "50"}],
    )
    assert result[0]["qual"] == "50"


def test_qual_missing_one_side():
    result = run_merge(
        [{"ref": "A", "alts": ["T"], "qual": "."}],
        [{"ref": "A", "alts": ["T"], "qual": "40"}],
    )
    assert result[0]["qual"] == "40"


def test_id_joined():
    result = run_merge(
        [{"ref": "A", "alts": ["T"], "id": "rs1"}],
        [{"ref": "A", "alts": ["T"], "id": "rs2"}],
    )
    assert result[0]["id"] == "rs1;rs2"


def test_filter_union():
    result = run_merge(
        [{"ref": "A", "alts": ["T"], "filt": "LowQual"}],
        [{"ref": "A", "alts": ["T"], "filt": "PASS"}],
    )
    assert result[0]["filt"] == "LowQual"


def test_ordering_preserved():
    result = run_merge(
        [{"ref": "A", "alts": ["T"]}, {"ref": "A", "alts": ["C"]}],
        [{"ref": "A", "alts": ["G"]}, {"ref": "A", "alts": ["C"]}],
    )
    assert [r["alts"] for r in result] == [["T"], ["G"], ["C"]]
