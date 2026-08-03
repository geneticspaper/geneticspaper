#!/usr/bin/env python3
"""Fetch founder-haplogroup metadata from GenBank/NCBI Nucleotide.

The core Mitotree workbook already includes many GenBank accessions and the
Behar/Costa studies, but this script adds a refreshable GenBank pull for the
four major Ashkenazi founder haplogroups. It uses NCBI E-utilities politely
(tool/email params, small request count, backoff on 429/5xx) and writes a CSV
that the R enrichment step merges and deduplicates with Mitotree S1.

Usage:
    python3 scripts/01_fetch_genbank_founders.py [email]

If no email is provided, a placeholder is sent; set one for production use.
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "external", "genbank_founder_records.csv")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FOUNDERS = ("K1a1b1a", "K1a9", "K2a2a", "N1b2")


def request(url: str, tries: int = 6) -> bytes:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mitomdoel/0.1"})
            with urllib.request.urlopen(req, timeout=60) as r:
                time.sleep(0.45)
                return r.read()
        except urllib.error.HTTPError as e:
            last = e
            if e.code not in (429, 500, 502, 503, 504):
                raise
            time.sleep(2 ** attempt)
        except Exception as e:  # network/transient
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"NCBI request failed after retries: {last}")


def esearch(term: str, email: str, retmax: int = 500) -> list[str]:
    params = {
        "db": "nuccore",
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
        "tool": "mitomdoel",
        "email": email,
    }
    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    data = json.loads(request(url).decode("utf-8"))
    return data["esearchresult"].get("idlist", [])


def esummary(ids: list[str], email: str) -> list[dict[str, str]]:
    if not ids:
        return []
    rows: list[dict[str, str]] = []
    for i in range(0, len(ids), 150):
        chunk = ids[i : i + 150]
        params = {
            "db": "nuccore",
            "id": ",".join(chunk),
            "retmode": "xml",
            "tool": "mitomdoel",
            "email": email,
        }
        url = f"{EUTILS}/esummary.fcgi?{urllib.parse.urlencode(params)}"
        root = ET.fromstring(request(url))
        for doc in root.findall(".//DocSum"):
            item = {child.attrib.get("Name", ""): (child.text or "") for child in doc.findall("Item")}
            rows.append(
                {
                    "uid": doc.findtext("Id", default=""),
                    "accession_version": item.get("AccessionVersion", ""),
                    "title": item.get("Title", ""),
                    "organism": item.get("Organism", ""),
                    "length": item.get("Length", ""),
                    "moltype": item.get("MolType", ""),
                    "subtype": item.get("SubType", ""),
                    "subname": item.get("SubName", ""),
                    "taxid": item.get("TaxId", ""),
                    "extra": item.get("Extra", ""),
                }
            )
    return rows


def main() -> None:
    email = sys.argv[1] if len(sys.argv) > 1 else "andrew@example.com"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    all_rows: list[dict[str, str]] = []
    seen = set()
    for hap in FOUNDERS:
        term = f'"{hap}"[All Fields] AND "Homo sapiens"[Organism] AND mitochondrion[All Fields]'
        print(f"Searching GenBank: {term}")
        ids = esearch(term, email=email)
        print(f"  {hap}: {len(ids)} IDs")
        for row in esummary(ids, email=email):
            acc = row["accession_version"] or row["uid"]
            key = (hap, acc)
            if key in seen:
                continue
            seen.add(key)
            row["query_haplogroup"] = hap
            row["query_term"] = term
            row["source"] = "GenBank_EUtilities"
            all_rows.append(row)
    fields = [
        "query_haplogroup", "source", "uid", "accession_version", "title",
        "organism", "length", "moltype", "subtype", "subname", "taxid",
        "extra", "query_term",
    ]
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)
    print(f"Wrote {len(all_rows)} GenBank rows -> {OUT}")


if __name__ == "__main__":
    main()
