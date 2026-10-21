#!/usr/bin/env python3
"""Fetch raw NCBI RunInfo metadata; never download reads or control samples.

The cached Nextflow task owns this network operation. Cohort parsing, scientific
eligibility, BioSample/run grouping, and scheduling are implemented in Groovy/DSL2.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def request(endpoint: str, parameters: dict[str, str]) -> bytes:
    payload = urllib.parse.urlencode(parameters).encode("ascii")
    for attempt in range(4):
        query = urllib.request.Request(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + endpoint,
            data=payload,
            headers={"User-Agent": "metagenomics-nextflow/3.0", "Content-Type": "application/x-www-form-urlencoded"},
        )
        try:
            with urllib.request.urlopen(query, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 3:
                raise
        except urllib.error.URLError:
            if attempt == 3:
                raise
        time.sleep(2**attempt)
    raise RuntimeError("NCBI metadata request failed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--email", default="")
    args = parser.parse_args()
    if not re.fullmatch(r"PRJ(?:NA|EB|DB)[0-9]+", args.project):
        parser.error("Invalid BioProject accession")
    common = {"db": "sra", "tool": "metagenomics_nextflow"}
    if args.email:
        common["email"] = args.email
    search = json.loads(request("esearch.fcgi", {
        **common, "term": f"{args.project}[BioProject]", "retmode": "json",
        "retmax": "0", "usehistory": "y",
    }))["esearchresult"]
    count = int(search["count"])
    if count < 1:
        raise ValueError(f"No SRA metadata found for {args.project}")
    raw = request("efetch.fcgi", {
        **common, "query_key": search["querykey"], "WebEnv": search["webenv"],
        "rettype": "runinfo", "retmode": "text", "retstart": "0", "retmax": str(count),
    })
    header = next(csv.reader(io.StringIO(raw.decode("utf-8-sig"))), [])
    if not raw.strip() or "Run" not in header:
        raise ValueError("NCBI did not return a RunInfo CSV")
    temporary = args.output.with_suffix(args.output.suffix + ".partial")
    temporary.write_bytes(raw)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
