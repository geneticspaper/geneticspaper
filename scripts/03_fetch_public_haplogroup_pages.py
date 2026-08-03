#!/usr/bin/env python3
"""Fetch public YFull / FTDNA haplogroup pages and extract lightweight metadata.

This is an enrichment layer, not a primary sampling source. It only uses public
pages, caches the raw HTML/text, and records parse status so missing or blocked
pages are explicit rather than silently hand-filled.
"""

from __future__ import annotations

import csv
import html
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "data" / "reference" / "public_haplogroup_page_targets.csv"
OUT_DIR = ROOT / "data" / "public_haplogroup_pages"
CACHE = OUT_DIR / "cache"
SUMMARY = OUT_DIR / "public_haplogroup_page_summary.csv"


def fetch(url: str, cache_path: Path, delay: float = 1.0) -> tuple[str, str]:
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8", errors="replace"), "cached"
    req = Request(
        url,
        headers={
            "User-Agent": "mitomdoel-research/1.0 (+public haplogroup metadata audit)",
            "Accept": "text/html,application/xhtml+xml,text/plain;q=0.9,*/*;q=0.8",
        },
    )
    try:
        with urlopen(req, timeout=30) as resp:
            text = resp.read().decode("utf-8", errors="replace")
        cache_path.write_text(text, encoding="utf-8")
        time.sleep(delay)
        return text, "fetched"
    except (HTTPError, URLError, TimeoutError) as exc:
        return "", f"error:{type(exc).__name__}:{exc}"


def strip_tags(text: str) -> str:
    text = re.sub(r"(?is)<script.*?</script>|<style.*?</style>", " ", text)
    text = re.sub(r"(?s)<[^>]+>", " ", text)
    text = html.unescape(text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def first(pattern: str, text: str, flags: int = re.I) -> str:
    m = re.search(pattern, text, flags)
    return m.group(1).strip() if m else ""


def all_join(pattern: str, text: str, flags: int = re.I, limit: int = 20) -> str:
    vals = []
    for m in re.finditer(pattern, text, flags):
        val = m.group(1).strip()
        if val and val not in vals:
            vals.append(val)
        if len(vals) >= limit:
            break
    return ";".join(vals)


def parse_yfull(raw: str) -> dict[str, str]:
    text = strip_tags(raw)
    ids = re.findall(r"\bid:([A-Z0-9_.-]+)", text)
    countries = re.findall(r"\[([A-Z]{2}(?:-[A-Z0-9]+)?)\]", text)
    return {
        "yfull_status": "ok" if text else "missing",
        "yfull_formed": first(r"formed\s+([0-9,]+\s+ybp)", text),
        "yfull_tmrca": first(r"TMRCA\s+([0-9,]+\s+ybp)", text),
        "yfull_public_ids_n": str(len(set(ids))) if ids else "",
        "yfull_country_tags": ";".join(sorted(set(countries)))[:300],
        "yfull_sample_ids": ";".join(ids[:20]),
    }


def parse_ftdna(raw: str) -> dict[str, str]:
    if not raw:
        return {
            "ftdna_status": "missing",
            "ftdna_formed_mean": "",
            "ftdna_formed_95ci": "",
            "ftdna_tmrca_mean": "",
            "ftdna_tmrca_95ci": "",
            "ftdna_children_n": "",
            "ftdna_placements_n": "",
            "ftdna_modern_country_counts": "",
            "ftdna_ancient_codes": "",
            "ftdna_ancient_studies": "",
            "ftdna_project_counts": "",
            "ftdna_variants_n": "",
        }
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        return {"ftdna_status": "unparsed_json"}

    def time_field(kind: str, field: str) -> str:
        value = data.get("time", {}).get(kind, {}).get(field)
        return "" if value is None else str(value)

    def ci(kind: str) -> str:
        oldest = time_field(kind, "oldest")
        youngest = time_field(kind, "youngest")
        return f"{oldest}..{youngest}" if oldest and youngest else ""

    countries = []
    for c in data.get("countries", [])[:12]:
        name = c.get("name") or ""
        count = c.get("count")
        total = c.get("total")
        if name and count is not None:
            denom = f"/{total}" if total else ""
            countries.append(f"{name}:{count}{denom}")

    ancient_codes = []
    ancient_studies = []
    for a in data.get("ancient", [])[:30]:
        code = a.get("code")
        if code:
            ancient_codes.append(code)
        for study in a.get("studies") or []:
            if study not in ancient_studies:
                ancient_studies.append(study)

    projects = []
    for p in data.get("projects", [])[:12]:
        name = p.get("ShortName") or p.get("Name") or ""
        members = p.get("HaplogroupPublicMembers")
        total = p.get("TotalPublicMembers")
        if name and members is not None and total is not None:
            projects.append(f"{name}:{members}/{total}")

    return {
        "ftdna_status": "ok",
        "ftdna_formed_mean": time_field("formed", "mean"),
        "ftdna_formed_95ci": ci("formed"),
        "ftdna_tmrca_mean": time_field("tmrca", "mean"),
        "ftdna_tmrca_95ci": ci("tmrca"),
        "ftdna_children_n": str(len(data.get("descendants", {}).get("children") or [])),
        "ftdna_placements_n": str(data.get("descendants", {}).get("placements") or ""),
        "ftdna_modern_country_counts": ";".join(countries),
        "ftdna_ancient_codes": ";".join(ancient_codes),
        "ftdna_ancient_studies": ";".join(ancient_studies),
        "ftdna_project_counts": ";".join(projects),
        "ftdna_variants_n": str(len(data.get("variants") or [])),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(TARGETS.open(newline="", encoding="utf-8")))
    out_rows = []
    for row in rows:
        key = row["lineage_key"]
        yfull_name = row["yfull_name"]
        ftdna_name = row["ftdna_name"]
        yfull_url = f"https://www.yfull.com/mtree/{yfull_name}/" if yfull_name else ""
        ftdna_story_url = f"https://discover.familytreedna.com/mtdna/{ftdna_name}/story" if ftdna_name else ""
        ftdna_url = f"https://discover.familytreedna.com/resources/mtdna/{ftdna_name}.json" if ftdna_name else ""

        y_raw, y_fetch = fetch(yfull_url, CACHE / f"yfull_{key}.html") if yfull_url else ("", "")
        f_raw, f_fetch = fetch(ftdna_url, CACHE / f"ftdna_{key}.json") if ftdna_url else ("", "")

        rec = {
            **row,
            "yfull_url": yfull_url,
            "ftdna_url": ftdna_url,
            "ftdna_story_url": ftdna_story_url,
            "yfull_fetch": y_fetch,
            "ftdna_fetch": f_fetch,
        }
        rec.update(parse_yfull(y_raw))
        rec.update(parse_ftdna(f_raw))
        out_rows.append(rec)

    fieldnames = list(out_rows[0].keys()) if out_rows else []
    with SUMMARY.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote {SUMMARY} ({len(out_rows)} rows)")


if __name__ == "__main__":
    main()
