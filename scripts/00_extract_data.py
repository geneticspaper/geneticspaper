#!/usr/bin/env python3
"""Stream-extract worksheets from the Mitotree supplementary .xlsx to CSV.

The workbook (media-1(1).xlsx) contains a 46 MB "S1 Public Samples" sheet that
makes ordinary loaders (openpyxl / readxl) slow or memory-hungry. This script
uses only the Python standard library (zipfile + regex) to stream the shared
string table and each worksheet row-by-row, writing plain CSV files that the R
pipeline consumes.

Usage:
    python3 scripts/00_extract_data.py [path/to/workbook.xlsx]

The workbook path may be passed as the first argument or configured via the
MITOTREE_XLSX key in a project-root `.env` file (see `.env.example`).
"""
from __future__ import annotations

import csv
import html
import os
import re
import sys
import zipfile

import _env

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT_XLSX = _env.get("MITOTREE_XLSX", "")
OUT_DIR = os.path.join(ROOT, "data", "derived")

# Sheets we want -> output CSV filename (by workbook sheet *name*).
WANTED = {
    "S1 Public Samples": "s1_public_samples.csv",
    "S10 Mitotree Structure": "s10_mitotree_structure.csv",
}

SI_RE = re.compile(rb"<si>(.*?)</si>", re.S)
T_RE = re.compile(rb"<t[^>]*>(.*?)</t>", re.S)
ROW_RE = re.compile(rb"<row[^>]*>(.*?)</row>", re.S)
# capture column letters, optional type attr, optional <v> value, optional inline <t>
CELL_RE = re.compile(
    rb'<c r="([A-Z]+)\d+"(?:[^>]*?t="(\w+)")?[^>]*?>'
    rb"(?:<v>(.*?)</v>|<is>(?:<t[^>]*>(.*?)</t>)?</is>)?",
    re.S,
)
SHEET_RE = re.compile(
    r'<sheet name="([^"]+)" sheetId="\d+"[^>]*r:id="(\w+)"'
)
REL_RE = re.compile(r'Id="(\w+)"[^>]*Target="(worksheets/sheet\d+\.xml)"')


def col_to_idx(col: str) -> int:
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - ord("A") + 1)
    return n - 1


def load_shared_strings(z: zipfile.ZipFile) -> list[str]:
    ss: list[str] = []
    with z.open("xl/sharedStrings.xml") as f:
        buf = b""
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            buf += chunk
            while True:
                m = SI_RE.search(buf)
                if not m:
                    break
                parts = T_RE.findall(m.group(1))
                ss.append(html.unescape(b"".join(parts).decode("utf-8", "replace")))
                buf = buf[m.end():]
    return ss


def sheet_file_map(z: zipfile.ZipFile) -> dict[str, str]:
    wb = z.read("xl/workbook.xml").decode("utf-8", "replace")
    rels = z.read("xl/_rels/workbook.xml.rels").decode("utf-8", "replace")
    relmap = dict(REL_RE.findall(rels))
    out: dict[str, str] = {}
    for name, rid in SHEET_RE.findall(wb):
        target = relmap.get(rid)
        if target:
            out[name] = "xl/" + target
    return out


def parse_row(rowxml: bytes, ss: list[str]) -> dict[int, str]:
    cells: dict[int, str] = {}
    for col, typ, v, inline in CELL_RE.findall(rowxml):
        col = col.decode()
        typ = typ.decode()
        idx = col_to_idx(col)
        if typ == "s" and v:
            cells[idx] = ss[int(v)]
        elif inline:
            cells[idx] = html.unescape(inline.decode("utf-8", "replace"))
        elif v:
            cells[idx] = html.unescape(v.decode("utf-8", "replace"))
        else:
            cells[idx] = ""
    return cells


def export_sheet(z: zipfile.ZipFile, path: str, ss: list[str], out_csv: str) -> int:
    n_written = 0
    header: list[str] | None = None
    ncols = 0
    with z.open(path) as f, open(out_csv, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        buf = b""
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            buf += chunk
            while True:
                m = ROW_RE.search(buf)
                if not m:
                    break
                cells = parse_row(m.group(1), ss)
                buf = buf[m.end():]
                if header is None:
                    ncols = (max(cells) + 1) if cells else 0
                    header = [cells.get(i, f"col{i}") for i in range(ncols)]
                    writer.writerow(header)
                    continue
                row = [cells.get(i, "") for i in range(ncols)]
                writer.writerow(row)
                n_written += 1
    return n_written


def main() -> None:
    xlsx = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_XLSX
    if not xlsx:
        sys.exit(
            "No workbook path given. Pass it as the first argument or set "
            "MITOTREE_XLSX in a project-root .env file (see .env.example)."
        )
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Reading workbook: {xlsx}")
    z = zipfile.ZipFile(xlsx)
    print("Loading shared strings ...")
    ss = load_shared_strings(z)
    print(f"  shared strings: {len(ss):,}")
    fmap = sheet_file_map(z)
    for name, out_name in WANTED.items():
        path = fmap.get(name)
        if not path:
            print(f"  WARNING: sheet {name!r} not found; skipping")
            continue
        out_csv = os.path.join(OUT_DIR, out_name)
        print(f"Exporting {name!r} -> {out_csv}")
        n = export_sheet(z, path, ss, out_csv)
        print(f"  wrote {n:,} data rows")


if __name__ == "__main__":
    main()
