# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///
"""Find sites with multiple records where at least one is multiallelic.

These are the sites where the bcftools anchor algorithm and connected-components
may disagree — disagreement requires 3+ records with at least one multiallelic.

Usage: uv run find_candidate_sites.py input.vcf[.gz]
"""

import subprocess
import sys
import tempfile
from itertools import groupby


def find_candidate_sites(vcf_path: str) -> list[tuple[str, str]]:
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%ALT\n", vcf_path],
        capture_output=True, text=True, check=True,
    )
    rows = (line.split("\t") for line in result.stdout.splitlines())
    candidates = []
    for (chrom, pos), group in groupby(rows, key=lambda r: (r[0], r[1])):
        records = list(group)
        if len(records) > 1 and any(len(r[2].split(",")) > 1 for r in records):
            candidates.append((chrom, pos))
    return candidates


def main(vcf_path: str) -> None:
    candidates = find_candidate_sites(vcf_path)
    if not candidates:
        return
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        for chrom, pos in candidates:
            f.write(f"{chrom}\t{pos}\n")
        sites_file = f.name
    subprocess.run(["bcftools", "view", "-R", sites_file, vcf_path], check=True)


if __name__ == "__main__":
    main(sys.argv[1])
