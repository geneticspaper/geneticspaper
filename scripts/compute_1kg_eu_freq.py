#!/usr/bin/env python3
"""Assign mtDNA haplogroups to the full 1000 Genomes Phase 3 European (EUR)
super-population and compute clade-specific non-Jewish frequency benchmarks for
the modeled lineages.

The EUR panel is all five 1000G European populations -- CEU, GBR (NW European),
FIN (Finnish), IBS (Iberian), TSI (Tuscan) -- n=503. Widening the panel beyond
the original CEU+GBR (n=190) roughly triples the non-Jewish sample size, tightens
the zero-frequency upper bounds for the founders, adds the Southern European
populations (Iberian/Tuscan) where Near-Eastern-affiliated control lineages
actually leak in, and lets several Near-Eastern control macro-backgrounds
(e.g. HV1b) resolve directly instead of climbing to an unrelated broader clade.

Uses callMom Phase 3 full mtDNA sequences (1000 Genomes Project Consortium, 2015)
+ Haplogrep3 (Phylotree 17 / rCRS). Related WES-based 1000G mtDNA work: Diroma et al.
2014 (BMC Genomics 15(Suppl 3):S2). Frequencies match the finest observed clade prefix
(longest match first), not macro-letter buckets.
"""

from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = ROOT / "data" / "external"
REFERENCE = ROOT / "data" / "reference"

PANEL = EXTERNAL / "1kg_phase3.panel"
FASTA = EXTERNAL / "chrMT_sequences.fasta.gz"
HAPLOGREP = EXTERNAL / "haplogrep3"
TARGETS = REFERENCE / "public_haplogroup_page_targets.csv"

WORK_FASTA = EXTERNAL / "eur_mt_clean.fasta"
HAPLOGREP_OUT = EXTERNAL / "eur_haplogroups.csv"

OUT_SAMPLES = REFERENCE / "1kg_eur_sample_haplogroups.csv"
OUT_FREQ = REFERENCE / "1kg_eur_clade_freq.csv"

# All five 1000G European populations (EUR super-population).
EUR_POPS = frozenset({"CEU", "GBR", "FIN", "IBS", "TSI"})
SOURCE = "1000 Genomes Phase 3 EUR (CEU+GBR+FIN+IBS+TSI, n={n})"
TREE = "phylotree-rcrs@17.3"
FOUNDERS = frozenset({"K1a1b1a", "K1a9", "K2a2a", "N1b2"})
MACRO_ROOT = {
    "V7a2c1b": "V",
    "U5a1f1a3": "U5",
    "L2a1l2a": "L2a",
    "M1a1b1c": "M1a",
    "M33c3": "M33",
    "N9a3a1b1": "N9a3",
    "A-a1b3a1": "A",
    "HV1b2": "HV1b",
    "R0a": "R0a",
    "U7a5": "U7",
    "U1b1a1": "U1",
    "H7": "H",
    "H6a1a1a": "H6",
    "J1c7a": "J1c",
    # Costa (2013) Suppl. Table S9/S7 published-origin calibration lineages.
    "T2": "T2",
    "W": "W",
    "I": "I",
    "U4": "U4",
    "H1": "H",
    "H3": "H",
    "H5": "H",
    "HV1a": "HV1",
}


def load_eur() -> dict[str, str]:
    samples: dict[str, str] = {}
    with PANEL.open() as fh:
        for line in fh:
            if line.startswith("sample"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in EUR_POPS:
                samples[parts[0]] = parts[1]
    return samples


def write_eur_fasta(eur: dict[str, str]) -> None:
    want = set(eur)
    with gzip.open(FASTA, "rt") as inf, WORK_FASTA.open("w") as outf:
        sid = None
        buf: list[str] = []
        for line in inf:
            if line.startswith(">"):
                if sid in want:
                    outf.write(f">{sid}\n")
                    outf.write("".join(buf).replace("-", ""))
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if sid in want:
            outf.write(f">{sid}\n")
            outf.write("".join(buf).replace("-", ""))


def run_haplogrep3() -> None:
    cmd = [
        str(HAPLOGREP),
        "classify",
        f"--in={WORK_FASTA}",
        f"--out={HAPLOGREP_OUT}",
        f"--tree={TREE}",
    ]
    subprocess.run(cmd, check=True, cwd=EXTERNAL)


def read_haplogrep(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = row.get("SampleID") or row.get("Sample") or ""
            hg = row.get("Haplogroup") or row.get("Haplogroup.main") or ""
            if sid and hg:
                rows.append((sid, hg))
    return rows


def load_target_aliases() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with TARGETS.open(newline="") as fh:
        for row in csv.DictReader(fh):
            key = row["lineage_key"]
            aliases = {key, row.get("ftdna_name", ""), row.get("yfull_name", "")}
            if key == "N1b2":
                aliases.add("N1b1b1")
            out[key] = {a for a in aliases if a}
    return out


def in_clade(hg: str, clade: str) -> bool:
    """Clade membership respecting mtDNA haplogroup boundaries.

    A subclade name continues with the opposite token type from the segment it
    extends (H1 -> H1a, L2a -> L2a1), so a character of the SAME token type
    immediately after the clade name denotes a *different* sibling clade, not a
    descendant: "H1" must not swallow H10/H11/H13, nor "H3" swallow H31, nor
    "H5" swallow H52; equally "H1a" must not swallow its siblings H1aa/H1ab
    (whose parent is H1+16189, not H1a), nor "H3a" swallow H3aa.

    PhyloTree also writes the node ancestral to X and Y as X'Y, so "H5'36" is
    the PARENT of H5 and H36, not a member of H5; a quote continuation is
    therefore rejected as well.
    """
    if hg == clade:
        return True
    if not hg.startswith(clade):
        return False
    nxt = hg[len(clade):len(clade) + 1]
    if nxt == "'":
        return False
    last = clade[-1:]
    if last.isdigit() and nxt.isdigit():
        return False
    if last.isalpha() and nxt.isalpha():
        return False
    return True


def count_clade(haplogroups: list[str], clade: str) -> int:
    return sum(1 for hg in haplogroups if in_clade(hg, clade))


def prefix_candidates(label: str, min_len: int = 2) -> list[str]:
    return [label[:i] for i in range(len(label), min_len - 1, -1)]


def best_clade_match(
    haplogroups: list[str],
    aliases: set[str],
    macro_root: str = "",
) -> tuple[str, int]:
    seen: set[str] = set()
    for alias in sorted(aliases, key=len, reverse=True):
        for pref in prefix_candidates(alias, min_len=2):
            if pref in seen:
                continue
            seen.add(pref)
            c = count_clade(haplogroups, pref)
            if c > 0:
                return pref, c
    if macro_root:
        for pref in prefix_candidates(macro_root, min_len=2):
            if pref in seen:
                continue
            seen.add(pref)
            c = count_clade(haplogroups, pref)
            if c > 0:
                return pref, c
        if len(macro_root) == 1:
            c = count_clade(haplogroups, macro_root)
            if c > 0:
                return macro_root, c
    return "", 0


def main() -> int:
    for path in (PANEL, FASTA, HAPLOGREP, TARGETS):
        if not path.exists():
            print(f"Missing required input: {path}", file=sys.stderr)
            return 1

    eur = load_eur()
    # Reuse a prior Haplogrep3 classification when present so frequency-only
    # refreshes (e.g. adding target lineages) do not require re-running the
    # Java classifier; delete eur_haplogroups.csv to force reclassification.
    if not HAPLOGREP_OUT.exists():
        write_eur_fasta(eur)
        run_haplogrep3()

    assigned = read_haplogrep(HAPLOGREP_OUT)
    pop_map = eur
    rows = [
        (sid, pop_map.get(sid, ""), hg)
        for sid, hg in assigned
        if sid in pop_map
    ]
    n = len(rows)
    if n != len(eur):
        missing = sorted(set(eur) - {sid for sid, _, _ in rows})
        print(f"Warning: classified {n}/{len(eur)} panel samples", file=sys.stderr)
        if missing[:5]:
            print(f"  missing e.g. {', '.join(missing[:5])}", file=sys.stderr)

    OUT_SAMPLES.parent.mkdir(parents=True, exist_ok=True)
    with OUT_SAMPLES.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "population", "haplogroup", "n_panel"])
        for sid, pop, hg in sorted(rows):
            w.writerow([sid, pop, hg, n])

    hap_counts: dict[str, int] = defaultdict(int)
    haps = [hg for _, _, hg in rows]
    for hg in haps:
        hap_counts[hg] += 1

    target_aliases = load_target_aliases()
    freq_rows = []
    for key, aliases in sorted(target_aliases.items()):
        if key in FOUNDERS:
            continue
        matched_clade, count = best_clade_match(
            haps, aliases, macro_root=MACRO_ROOT.get(key, "")
        )
        freq_rows.append({
            "haplotype": matched_clade or key,
            "lineage_key": key,
            "count": count,
            "frequency": count / n if n else 0.0,
            "n_nonjews": n,
            "source": SOURCE.format(n=n),
        })

    with OUT_FREQ.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["haplotype", "lineage_key", "count", "frequency", "n_nonjews", "source"],
        )
        w.writeheader()
        w.writerows(freq_rows)

    print(f"Assigned {n} EUR samples -> {OUT_SAMPLES}")
    print(f"Wrote {len(freq_rows)} clade frequencies -> {OUT_FREQ}")
    top = sorted(hap_counts.items(), key=lambda x: -x[1])[:12]
    print("Top assigned haplogroups:", ", ".join(f"{h}={c}" for h, c in top))
    for key in ("V7a2c1b", "U5a1f1a3", "H7", "HV1b2", "J1c7a"):
        hit = next((r for r in freq_rows if r["lineage_key"] == key), None)
        if hit:
            print(
                f"  {key}: {hit['count']}/{n} ({100 * hit['frequency']:.2f}%) "
                f"via clade {hit['haplotype'] or '(none)'}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
