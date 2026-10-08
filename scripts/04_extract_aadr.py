#!/usr/bin/env python3
"""Extract ancient mtDNA calls from the Allen Ancient DNA Resource (AADR) v66
annotation file into a compact CSV for the ancient-DNA falsification analysis.

The AADR (Mallick et al. 2024) is the large public compilation of ancient and
present-day genotypes; its `.anno` metadata file records, per individual, a
manually/published mtDNA haplogroup call, a mean date in years before 1950,
geographic coordinates, political entity, and a free-text group label. This
script keeps only the fields the origin analysis needs and writes them to
`data/external/aadr_mtdna.csv`.

The annotation file is large and is not committed; download it from the AADR
Harvard Dataverse (doi:10.7910/DVN/FFIDCW) and either place it at the default
path below or set AADR_ANNO in a project-root `.env` file (see `.env.example`).
"""
from __future__ import annotations

import csv
import os
import sys

import _env

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT_ANNO = _env.get(
    "AADR_ANNO", os.path.join(ROOT, "data", "external", "aadr_v66_1240K.anno")
)
OUT = os.path.join(ROOT, "data", "external", "aadr_mtdna.csv")

# 1-based column indices in the AADR v66 `.anno` file.
COL_ID = 1
COL_DATE_BP = 11        # mean date in years before 1950 CE (0 = present-day)
COL_GROUP = 15          # Group ID (free-text label, e.g. Germany_Medieval_Jewish)
COL_COUNTRY = 17        # Political Entity
COL_LAT = 18
COL_LON = 19
COL_MT = 39             # mtDNA haplogroup if >2x or published

# A few AADR political-entity spellings the R region classifier does not key on.
COUNTRY_ALIAS = {
    "Czechia": "Czech Republic",
    "Türkiye": "Turkey",
    "Turkiye": "Turkey",
    "Russia": "Russian Federation",
    "Palestine": "Palestinian Territory",
}

JEWISH_TOKENS = ("jewish", "jew_", "ashkenaz")


def parse_float(x: str):
    try:
        v = float(x)
        return v
    except (ValueError, TypeError):
        return ""


def clean_mt(hg: str) -> str:
    hg = (hg or "").strip()
    if not hg:
        return ""
    # Some calls carry a '+' extra-mutation suffix or trailing coverage note;
    # keep the haplogroup token itself (first whitespace-delimited field).
    base = hg.split()[0]
    # AADR writes missing/undetermined calls as '..', '.', or 'n/a (<2x)' -- note
    # the coverage note, so the test must come AFTER taking the first token, or
    # 'n/a' survives as if it were a haplogroup name. 'rCRS' is a reference
    # placeholder, not a call.
    if base.lower() in ("", "..", ".", "n/a", "na", "rcrs"):
        return ""
    return base


def main() -> int:
    anno = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ANNO
    if not anno or not os.path.exists(anno):
        sys.exit(
            f"AADR annotation file not found ({anno!r}). Download it from the AADR "
            "Harvard Dataverse (doi:10.7910/DVN/FFIDCW) and set AADR_ANNO in .env "
            "or pass the path as the first argument."
        )

    n_in = n_out = 0
    with open(anno, encoding="utf-8", errors="replace") as fh, open(
        OUT, "w", newline="", encoding="utf-8"
    ) as out:
        w = csv.writer(out)
        w.writerow(
            ["id", "mt_haplogroup", "date_bp", "year", "country", "group_id",
             "is_jewish", "lat", "lon"]
        )
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            n_in += 1
            if len(row) <= COL_MT - 1:
                continue
            mt = clean_mt(row[COL_MT - 1])
            if not mt:
                continue
            date_bp = parse_float(row[COL_DATE_BP - 1])
            year = "" if date_bp == "" else int(round(1950 - date_bp))
            country = row[COL_COUNTRY - 1].strip()
            country = COUNTRY_ALIAS.get(country, country)
            group = row[COL_GROUP - 1].strip()
            is_jewish = int(any(t in group.lower() for t in JEWISH_TOKENS))
            w.writerow([
                row[COL_ID - 1].strip(), mt,
                "" if date_bp == "" else date_bp, year,
                country, group, is_jewish,
                parse_float(row[COL_LAT - 1]), parse_float(row[COL_LON - 1]),
            ])
            n_out += 1

    print(f"Read {n_in} AADR rows; wrote {n_out} with mtDNA calls -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
