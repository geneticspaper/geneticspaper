#!/usr/bin/env python3
"""Extract useful Genes 2026 supplement tables for the mtDNA report.

Inputs are the MDPI Genes 2026 Tàrrega/Catalonia medieval Jewish aDNA paper:
  - genes-17-00358.xml
  - genes-17-00358-s001/Supplementary_material_TablesS1-S15.xlsx

Outputs:
  - data/external/genes2026_samples.csv
  - data/external/genes2026_mtdna.csv
  - data/external/genes2026_qpadm.csv

Only the standard library is used so this remains reproducible without pandas.

The workbook path may be passed as the first argument or configured via the
GENES2026_XLSX key in a project-root `.env` file (see `.env.example`).
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
DEFAULT_XLSX = _env.get("GENES2026_XLSX", "")
OUT_DIR = os.path.join(ROOT, "data", "external")

WANTED = {
    "Table S1 - Samples": "genes2026_samples.csv",
    "Table S6 - mtDNA haplogroups": "genes2026_mtdna.csv",
    "Table S15 - qpAdm": "genes2026_qpadm.csv",
}


def col_to_idx(col: str) -> int:
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - ord("A") + 1)
    return n - 1


def shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    data = z.read("xl/sharedStrings.xml")
    out = []
    for si in re.findall(rb"<si>(.*?)</si>", data, re.S):
        parts = re.findall(rb"<t[^>]*>(.*?)</t>", si, re.S)
        out.append(html.unescape(b"".join(parts).decode("utf-8", "replace")))
    return out


def sheet_map(z: zipfile.ZipFile) -> dict[str, str]:
    wb = z.read("xl/workbook.xml").decode("utf-8", "replace")
    rels = z.read("xl/_rels/workbook.xml.rels").decode("utf-8", "replace")
    relmap = dict(re.findall(r'Id="(\w+)"[^>]*Target="(worksheets/sheet\d+\.xml)"', rels))
    out = {}
    for name, _sid, rid in re.findall(r'<sheet name="([^"]+)" sheetId="(\d+)" r:id="(\w+)"', wb):
        if rid in relmap:
            out[name] = "xl/" + relmap[rid]
    return out


CELL_RE = re.compile(
    rb'<c r="([A-Z]+)\d+"(?:[^>]*?t="(\w+)")?[^>]*?>'
    rb"(?:<v>(.*?)</v>|<is>(.*?)</is>)?",
    re.S,
)
TEXT_RE = re.compile(rb"<t[^>]*>(.*?)</t>", re.S)


def rows_from_sheet(z: zipfile.ZipFile, path: str, ss: list[str]) -> list[list[str]]:
    xml = z.read(path)
    rows = []
    for rowxml in re.findall(rb"<row[^>]*>(.*?)</row>", xml, re.S):
        cells = {}
        for col, typ, val, inline in CELL_RE.findall(rowxml):
            idx = col_to_idx(col.decode())
            typ = typ.decode()
            value = val.decode() if val else ""
            if typ == "s" and value:
                text = ss[int(value)]
            elif inline:
                text = html.unescape(b"".join(TEXT_RE.findall(inline)).decode("utf-8", "replace"))
            else:
                text = html.unescape(value)
            cells[idx] = re.sub(r"\s+", " ", text).strip()
        if cells:
            rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
    return rows


def find_header(rows: list[list[str]], must_contain: str) -> int:
    for i, row in enumerate(rows):
        if any(must_contain.lower() == cell.lower() for cell in row):
            return i
    raise ValueError(f"Header containing {must_contain!r} not found")


def export_sheet(rows: list[list[str]], header_idx: int, out_path: str) -> int:
    header = rows[header_idx]
    # Drop trailing empty columns.
    while header and not header[-1]:
        header = header[:-1]
    ncols = len(header)
    data = []
    for row in rows[header_idx + 1 :]:
        row = (row + [""] * ncols)[:ncols]
        if not any(row):
            continue
        # Stop at fully blank/table-note areas.
        data.append(row)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)
    return len(data)


def main() -> None:
    xlsx = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_XLSX
    if not xlsx:
        sys.exit(
            "No workbook path given. Pass it as the first argument or set "
            "GENES2026_XLSX in a project-root .env file (see .env.example)."
        )
    os.makedirs(OUT_DIR, exist_ok=True)
    z = zipfile.ZipFile(xlsx)
    ss = shared_strings(z)
    smap = sheet_map(z)
    headers = {
        "Table S1 - Samples": "Genetic code",
        "Table S6 - mtDNA haplogroups": "Sample",
        "Table S15 - qpAdm": "Target",
    }
    for sheet, out_name in WANTED.items():
        rows = rows_from_sheet(z, smap[sheet], ss)
        idx = find_header(rows, headers[sheet])
        out = os.path.join(OUT_DIR, out_name)
        n = export_sheet(rows, idx, out)
        print(f"{sheet}: wrote {n} rows -> {out}")


if __name__ == "__main__":
    main()
