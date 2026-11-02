#!/usr/bin/env python3
"""
generate_pipeline_summaries.py — Modular pipeline summaries generator.

Produces three conceptually separated summaries:
  1. taxonomy_summary.tsv: Read-based taxonomic results from Kraken2 and Bracken.
  2. mag_summary.tsv: MAG-based results from CheckM2, filtering, GTDB-Tk, CoverM.
  3. run_overview.tsv: Technical run overview, stages, inputs, configuration, inventories.

Exit codes:
  0: Summary created successfully.
  1: Source present but invalid or inconsistent.
  2: Usage or configuration error.
  3: Mandatory source absent.
  4: Summary not-applicable.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

SCHEMA_VERSION = "1"

STATUS_ORDER = [
    "quality-reported",
    "retained",
    "excluded",
    "classified",
    "unclassified",
    "coverm-reported",
    "coverm-zero",
    "checkm2-only",
]

ALLOWED_COVERM_METRICS = {
    "relative_abundance",
    "covered_fraction",
    "mean",
    "mean_coverage",
    "variance",
    "trimmed_mean",
    "covered_bases",
    "count",
    "reads_per_base",
    "rpkm",
    "tpm",
}


def sanitize_tsv(val):
    if val is None:
        return "NA"
    s = str(val).strip()
    if not s:
        return "NA"
    return s.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def compute_sha256(path):
    if not path or not os.path.isfile(path):
        return "unknown"
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            while chunk := f.read(65536):
                h.update(chunk)
        return h.hexdigest().lower()
    except OSError:
        return "unknown"


def parse_sample_manifest(path):
    samples = []
    if not path or not os.path.isfile(path):
        return samples
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for i, row in enumerate(reader):
            if not row or not any(row):
                continue
            if i == 0 and ("sample" in row[0].lower() or "sample_id" in row[0].lower()):
                continue
            sample = row[0].strip()
            if sample and sample not in samples:
                samples.append(sample)
    return samples


def atomic_publish_file(source_path, target_path):
    target_dir = os.path.dirname(os.path.abspath(target_path))
    os.makedirs(target_dir, exist_ok=True)
    temp_target = os.path.join(
        target_dir, f".atomic.{os.getpid()}.{os.path.basename(target_path)}"
    )
    shutil.copy2(source_path, temp_target)
    os.replace(temp_target, target_path)


# ------------------------------------------------------------------------------
# Taxonomy Summary Subcommand
# ------------------------------------------------------------------------------
def run_taxonomy(args):
    if not args.kraken2_enabled:
        # Not applicable when Kraken2 is disabled
        if args.dry_run:
            print("Taxonomy summary : not-applicable")
        sys.exit(4)

    samples = parse_sample_manifest(args.sample_manifest)
    if not samples and not args.dry_run:
        print("[ERROR] No valid samples found in sample manifest.", file=sys.stderr)
        sys.exit(3)

    kraken_dir = args.kraken2_dir
    if not kraken_dir or not os.path.isdir(kraken_dir):
        if args.dry_run:
            print("Taxonomy summary : would run (Kraken2 directory deferred)")
            sys.exit(0)
        print(f"[ERROR] Kraken2 directory missing: {kraken_dir}", file=sys.stderr)
        sys.exit(3)

    records = []
    seen_keys = set()
    source_files = []
    source_fps = {}

    for sample in samples:
        # 1. Kraken2 report
        k2_candidates = [
            os.path.join(kraken_dir, f"kraken_report_{sample}.txt"),
            os.path.join(kraken_dir, f"{sample}.kraken2.report.txt"),
            os.path.join(kraken_dir, f"{sample}.kraken_report.txt"),
            os.path.join(kraken_dir, f"kraken_{sample}.report.txt"),
            os.path.join(kraken_dir, f"{sample}.report.txt"),
        ]
        k2_report = None
        for cand in k2_candidates:
            if os.path.isfile(cand):
                k2_report = cand
                break

        if not k2_report:
            if args.dry_run:
                continue
            print(f"[ERROR] Missing Kraken2 report for sample {sample}", file=sys.stderr)
            sys.exit(3)

        if os.path.getsize(k2_report) == 0:
            if not args.dry_run:
                print(f"[ERROR] Empty Kraken2 report for sample {sample}: {k2_report}", file=sys.stderr)
                sys.exit(1)

        source_files.append(k2_report)
        source_fps[k2_report] = compute_sha256(k2_report)

        with open(k2_report, "r", encoding="utf-8") as f:
            for line_no, line in enumerate(f, 1):
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 6:
                    print(f"[ERROR] Inconsistent Kraken2 line {line_no} in {k2_report}", file=sys.stderr)
                    sys.exit(1)
                pct_str, clade_str, direct_str, rank_code, taxid_str, raw_name = parts[:6]
                try:
                    pct = float(pct_str)
                    direct_reads = int(direct_str)
                    taxid = int(taxid_str)
                except ValueError:
                    print(f"[ERROR] Invalid numeric value in Kraken2 report {k2_report}:{line_no}", file=sys.stderr)
                    sys.exit(1)

                fraction = pct / 100.0
                if fraction < 0.0 or fraction > 1.0:
                    print(f"[ERROR] Kraken2 fraction out of range (0..1) in {k2_report}:{line_no}: {fraction}", file=sys.stderr)
                    sys.exit(1)

                taxon_name = raw_name.lstrip()
                taxonomic_level = rank_code.strip()

                if taxid == 0:
                    status = "unclassified"
                elif direct_reads == 0:
                    status = "zero-count"
                else:
                    status = "classified"

                key = (sample, "kraken2", taxonomic_level, taxid)
                if key in seen_keys:
                    print(f"[ERROR] Duplicate key in Kraken2 taxonomy summary: {key}", file=sys.stderr)
                    sys.exit(1)
                seen_keys.add(key)

                records.append({
                    "sample_id": sample,
                    "source": "kraken2",
                    "taxonomic_level": taxonomic_level,
                    "taxon_name": sanitize_tsv(taxon_name),
                    "taxonomy_id": taxid,
                    "assigned_reads": direct_reads,
                    "estimated_reads": "NA",
                    "fraction_total_reads": f"{fraction:.6g}",
                    "status": status,
                    "sort_taxid": taxid,
                })

        # 2. Bracken report (if enabled)
        if args.bracken_enabled:
            bracken_candidates = [
                os.path.join(kraken_dir, f"bracken_{sample}.txt"),
                os.path.join(kraken_dir, f"bracken_{sample}.tsv"),
                os.path.join(kraken_dir, f"{sample}.bracken.txt"),
                os.path.join(kraken_dir, f"{sample}.bracken.tsv"),
                os.path.join(kraken_dir, f"{sample}.bracken.{args.bracken_level}.tsv"),
                os.path.join(kraken_dir, f"{sample}.bracken.{args.bracken_level}.txt"),
            ]
            bracken_report = None
            for cand in bracken_candidates:
                if os.path.isfile(cand):
                    bracken_report = cand
                    break

            if not bracken_report:
                if args.dry_run:
                    continue
                print(f"[ERROR] Missing Bracken report for sample {sample}", file=sys.stderr)
                sys.exit(3)

            if os.path.getsize(bracken_report) == 0:
                if not args.dry_run:
                    print(f"[ERROR] Empty Bracken report for sample {sample}: {bracken_report}", file=sys.stderr)
                    sys.exit(1)

            source_files.append(bracken_report)
            source_fps[bracken_report] = compute_sha256(bracken_report)

            with open(bracken_report, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                fn_map = {fn.lower().strip(): fn for fn in (reader.fieldnames or [])}
                name_col = fn_map.get("name")
                taxid_col = fn_map.get("taxonomy_id") or fn_map.get("tax_id") or fn_map.get("taxid")
                level_col = fn_map.get("taxonomy_lvl") or fn_map.get("tax_code") or fn_map.get("level")
                assigned_col = fn_map.get("kraken_assigned_reads") or fn_map.get("assigned_reads")
                est_col = fn_map.get("new_est_reads") or fn_map.get("estimated_reads")
                frac_col = fn_map.get("fraction_total_reads") or fn_map.get("fraction")

                if not (name_col and taxid_col and level_col and assigned_col and est_col and frac_col):
                    print(f"[ERROR] Incomplete Bracken header in {bracken_report}", file=sys.stderr)
                    sys.exit(1)

                for line_no, row in enumerate(reader, 2):
                    try:
                        taxid = int(row[taxid_col].strip())
                        assigned = int(row[assigned_col].strip())
                        estimated = int(row[est_col].strip())
                        frac = float(row[frac_col].strip())
                    except (ValueError, KeyError):
                        print(f"[ERROR] Invalid numeric field in Bracken report {bracken_report}:{line_no}", file=sys.stderr)
                        sys.exit(1)

                    if frac < 0.0 or frac > 1.0:
                        print(f"[ERROR] Bracken fraction out of range (0..1) in {bracken_report}:{line_no}: {frac}", file=sys.stderr)
                        sys.exit(1)

                    taxon_name = row[name_col].strip()
                    taxonomic_level = row[level_col].strip()

                    status = "zero-count" if estimated == 0 else "estimated"

                    key = (sample, "bracken", taxonomic_level, taxid)
                    if key in seen_keys:
                        print(f"[ERROR] Duplicate key in Bracken taxonomy summary: {key}", file=sys.stderr)
                        sys.exit(1)
                    seen_keys.add(key)

                    records.append({
                        "sample_id": sample,
                        "source": "bracken",
                        "taxonomic_level": taxonomic_level,
                        "taxon_name": sanitize_tsv(taxon_name),
                        "taxonomy_id": taxid,
                        "assigned_reads": assigned,
                        "estimated_reads": estimated,
                        "fraction_total_reads": f"{frac:.6g}",
                        "status": status,
                        "sort_taxid": taxid,
                    })

    if args.dry_run:
        print(f"Taxonomy summary : would run ({len(records)} rows previewed)")
        sys.exit(0)

    # Sort deterministically: sample_id, source, taxonomic_level, taxonomy_id, taxon_name
    records.sort(key=lambda r: (r["sample_id"], r["source"], r["taxonomic_level"], r["sort_taxid"], r["taxon_name"]))

    # Write atomically
    out_dir = args.outdir
    os.makedirs(out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out_dir, prefix=f".build.taxonomy.{os.getpid()}_") as bdir:
        tmp_tsv = os.path.join(bdir, "taxonomy_summary.tsv")
        headers = [
            "sample_id",
            "source",
            "taxonomic_level",
            "taxon_name",
            "taxonomy_id",
            "assigned_reads",
            "estimated_reads",
            "fraction_total_reads",
            "status",
        ]
        with open(tmp_tsv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for r in records:
                writer.writerow(r)

        summary_fp = compute_sha256(tmp_tsv)

        # Write metadata JSON
        meta = {
            "schema_version": SCHEMA_VERSION,
            "summary_type": "taxonomy",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "pipeline_version": args.pipeline_version,
            "row_count": len(records),
            "source_files": source_files,
            "source_fingerprints": source_fps,
            "summary_fingerprint": summary_fp,
            "status": "available",
            "active_sample_manifest": args.sample_manifest,
            "sample_manifest_fingerprint": compute_sha256(args.sample_manifest),
            "kraken2_used": bool(args.kraken2_enabled),
            "bracken_used": bool(args.bracken_enabled),
            "bracken_level": args.bracken_level if args.bracken_enabled else "NA",
        }
        tmp_meta = os.path.join(bdir, "taxonomy_summary.meta.json")
        with open(tmp_meta, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)

        # Atomic publish
        final_tsv = os.path.join(out_dir, "taxonomy_summary.tsv")
        final_meta = os.path.join(out_dir, "taxonomy_summary.meta.json")
        atomic_publish_file(tmp_tsv, final_tsv)
        atomic_publish_file(tmp_meta, final_meta)

    print(f"[Taxonomy Summary] Published {len(records)} row(s) to {final_tsv}")
    sys.exit(0)


# ------------------------------------------------------------------------------
# MAG Summary Subcommand
# ------------------------------------------------------------------------------
def build_mag_status(tokens):
    # Order tokens in canonical fixed sequence
    active = [t for t in STATUS_ORDER if t in tokens]
    return ";".join(active) if active else "unknown"


def parse_gtdb_lineage(lineage_str):
    res = {
        "domain": "NA",
        "class": "NA",
        "order": "NA",
        "family": "NA",
        "genus": "NA",
        "species": "NA",
    }
    if not lineage_str or lineage_str in ("Unclassified", "NA", "None"):
        return res
    parts = lineage_str.split(";")
    for p in parts:
        p = p.strip()
        if p.startswith("d__"):
            val = p[3:].strip()
            res["domain"] = val if val else "NA"
        elif p.startswith("c__"):
            val = p[3:].strip()
            res["class"] = val if val else "NA"
        elif p.startswith("o__"):
            val = p[3:].strip()
            res["order"] = val if val else "NA"
        elif p.startswith("f__"):
            val = p[3:].strip()
            res["family"] = val if val else "NA"
        elif p.startswith("g__"):
            val = p[3:].strip()
            res["genus"] = val if val else "NA"
        elif p.startswith("s__"):
            val = p[3:].strip()
            res["species"] = val if val else "NA"
    return res


def run_mag(args):
    has_mag_source = (
        args.checkm2_enabled
        or args.filtering_enabled
        or args.gtdbtk_enabled
        or args.coverm_enabled
    )
    if not has_mag_source:
        if args.dry_run:
            print("MAG summary      : not-applicable")
        sys.exit(4)

    source_files = []
    source_fps = {}

    # CheckM2 report
    checkm2_data = {}  # mag_id -> (completeness, contamination)
    if args.checkm2_enabled:
        report_path = args.checkm2_report
        if not report_path or not os.path.isfile(report_path):
            if args.dry_run:
                print("MAG summary      : would run (CheckM2 deferred)")
                sys.exit(0)
            print(f"[ERROR] Missing CheckM2 quality report: {report_path}", file=sys.stderr)
            sys.exit(3)
        if os.path.getsize(report_path) == 0:
            print(f"[ERROR] Empty CheckM2 quality report: {report_path}", file=sys.stderr)
            sys.exit(1)

        source_files.append(report_path)
        source_fps[report_path] = compute_sha256(report_path)

        with open(report_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if not reader.fieldnames:
                print("[ERROR] CheckM2 report has no header.", file=sys.stderr)
                sys.exit(1)
            name_col = "Name" if "Name" in reader.fieldnames else ("Bin_ID" if "Bin_ID" in reader.fieldnames else reader.fieldnames[0])
            for line_no, row in enumerate(reader, 2):
                mag_id = row.get(name_col, "").strip()
                if not mag_id:
                    continue
                try:
                    comp = float(row.get("Completeness", "0").strip())
                    cont = float(row.get("Contamination", "0").strip())
                except ValueError:
                    print(f"[ERROR] Invalid numeric quality value in CheckM2 report:{line_no}", file=sys.stderr)
                    sys.exit(1)
                checkm2_data[mag_id] = (comp, cont)

    # Filtering manifest
    filter_data = {}  # mag_id -> (retained: bool, reason: str)
    if args.filtering_enabled:
        filt_path = args.filter_manifest
        if not filt_path or not os.path.isfile(filt_path):
            if args.dry_run:
                print("MAG summary      : would run (filtering manifest deferred)")
                sys.exit(0)
            print(f"[ERROR] Missing MAG filtering manifest: {filt_path}", file=sys.stderr)
            sys.exit(3)

        source_files.append(filt_path)
        source_fps[filt_path] = compute_sha256(filt_path)

        with open(filt_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            req = {"mag_id", "retained", "reason"}
            if not reader.fieldnames or not req.issubset(set(reader.fieldnames)):
                print("[ERROR] MAG filtering manifest missing required columns (mag_id, retained, reason).", file=sys.stderr)
                sys.exit(1)
            for row in reader:
                mid = row["mag_id"].strip()
                ret = row["retained"].strip().lower() == "true"
                reason = row["reason"].strip()
                filter_data[mid] = (ret, reason)

    # Base MAG IDs
    base_mags = list(checkm2_data.keys())
    if not base_mags and filter_data:
        base_mags = list(filter_data.keys())
    base_mags.sort()

    if not base_mags and not args.dry_run:
        print("[ERROR] No MAG IDs found in quality report or filtering manifest.", file=sys.stderr)
        sys.exit(1)

    # GTDB-Tk
    gtdb_data = {}  # mag_id -> dict(...)
    if args.gtdbtk_enabled:
        gtdb_files = []
        if args.gtdbtk_dir and os.path.isdir(args.gtdbtk_dir):
            for fname in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
                fp = os.path.join(args.gtdbtk_dir, fname)
                if os.path.isfile(fp) and os.path.getsize(fp) > 0:
                    gtdb_files.append(fp)
        if not gtdb_files and not args.dry_run:
            print(f"[ERROR] No valid GTDB-Tk summary TSV files found in {args.gtdbtk_dir}", file=sys.stderr)
            sys.exit(3)

        for gpath in gtdb_files:
            source_files.append(gpath)
            source_fps[gpath] = compute_sha256(gpath)
            with open(gpath, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if not reader.fieldnames or "user_genome" not in reader.fieldnames:
                    continue
                for line_no, row in enumerate(reader, 2):
                    mid = row.get("user_genome", "").strip()
                    if not mid:
                        continue
                    mid_norm = re.sub(r"\.(fa|fasta|fna)$", "", mid, flags=re.IGNORECASE)
                    clf = row.get("classification", "").strip() or "Unclassified"
                    method = row.get("classification_method", "").strip() or "NA"
                    closest_ref = row.get("closest_placement_reference", row.get("closest_reference", "NA")).strip() or "NA"
                    if closest_ref in ("N/A", "None", ""):
                        closest_ref = "NA"

                    ani_str = row.get("closest_placement_ani", row.get("closest_reference_ani", row.get("ani", "NA"))).strip()
                    af_str = row.get("closest_placement_af", row.get("closest_reference_af", row.get("af", "NA"))).strip()

                    ani_val = "NA"
                    if ani_str and ani_str not in ("NA", "N/A", "None"):
                        try:
                            ani_f = float(ani_str)
                            if ani_f < 0.0 or ani_f > 100.0 or math.isnan(ani_f):
                                print(f"[ERROR] GTDB-Tk ANI out of range (0..100) in {gpath}:{line_no}: {ani_f}", file=sys.stderr)
                                sys.exit(1)
                            ani_val = f"{ani_f:.4f}"
                        except ValueError:
                            print(f"[ERROR] Invalid GTDB-Tk ANI in {gpath}:{line_no}: {ani_str}", file=sys.stderr)
                            sys.exit(1)

                    af_val = "NA"
                    if af_str and af_str not in ("NA", "N/A", "None"):
                        try:
                            af_f = float(af_str)
                            if af_f < 0.0 or af_f > 1.0 or math.isnan(af_f):
                                print(f"[ERROR] GTDB-Tk AF out of range (0..1) in {gpath}:{line_no}: {af_f}", file=sys.stderr)
                                sys.exit(1)
                            af_val = f"{af_f:.4f}"
                        except ValueError:
                            print(f"[ERROR] Invalid GTDB-Tk AF in {gpath}:{line_no}: {af_str}", file=sys.stderr)
                            sys.exit(1)

                    gtdb_data[mid_norm] = {
                        "classification": clf,
                        "method": method,
                        "closest_reference": closest_ref,
                        "ani": ani_val,
                        "af": af_val,
                        "lineage": parse_gtdb_lineage(clf),
                    }

    # CoverM
    samples = parse_sample_manifest(args.sample_manifest)
    coverm_data = {}  # mag_id -> {sample_id -> {metric -> value}}
    coverm_files = []
    if args.coverm_enabled:
        cov_dir = args.coverm_dir
        if not cov_dir or not os.path.isdir(cov_dir):
            if args.dry_run:
                print("MAG summary      : would run (CoverM directory deferred)")
                sys.exit(0)
            print(f"[ERROR] CoverM directory missing: {cov_dir}", file=sys.stderr)
            sys.exit(3)

        per_sample_found = True
        temp_coverm_files = []
        if samples:
            for s in samples:
                candidates = [
                    os.path.join(cov_dir, f"abundance_{s}.tsv"),
                    os.path.join(cov_dir, f"{s}.abundance.tsv"),
                    os.path.join(cov_dir, f"coverm_{s}.tsv"),
                    os.path.join(cov_dir, f"{s}_abundance.tsv"),
                ]
                found = None
                for cand in candidates:
                    if os.path.isfile(cand):
                        found = cand
                        break
                if found:
                    temp_coverm_files.append((s, found))
                else:
                    per_sample_found = False
                    break
        else:
            per_sample_found = False

        if per_sample_found and temp_coverm_files:
            for s, cpath in temp_coverm_files:
                coverm_files.append((s, cpath))
                source_files.append(cpath)
                source_fps[cpath] = compute_sha256(cpath)
                with open(cpath, "r", encoding="utf-8") as f:
                    reader = csv.reader(f, delimiter="\t")
                    header = next(reader, None)
                    if not header or len(header) < 2:
                        continue
                    col_metrics = []
                    for col_name in header[1:]:
                        norm_col = col_name.lower().replace(s.lower(), "").strip()
                        norm_col = re.sub(r"[\s\-_]+", "_", norm_col).strip("_")
                        norm_col = norm_col.replace("(%)", "").strip("_")
                        if "relative_abundance" in norm_col:
                            canonical = "relative_abundance"
                        elif "covered_fraction" in norm_col:
                            canonical = "covered_fraction"
                        elif "mean_coverage" in norm_col:
                            canonical = "mean_coverage"
                        elif "trimmed_mean" in norm_col:
                            canonical = "trimmed_mean"
                        elif "covered_bases" in norm_col:
                            canonical = "covered_bases"
                        elif "reads_per_base" in norm_col:
                            canonical = "reads_per_base"
                        elif norm_col in ALLOWED_COVERM_METRICS:
                            canonical = norm_col
                        else:
                            canonical = norm_col
                        col_metrics.append(canonical)

                    for row in reader:
                        if not row or not row[0].strip():
                            continue
                        mid = row[0].strip()
                        mid_norm = re.sub(r"\.(fa|fasta|fna)$", "", mid, flags=re.IGNORECASE)
                        if mid_norm not in coverm_data:
                            coverm_data[mid_norm] = {}
                        if s not in coverm_data[mid_norm]:
                            coverm_data[mid_norm][s] = {}
                        for metric, val_str in zip(col_metrics, row[1:]):
                            try:
                                val = float(val_str.strip())
                                if math.isnan(val) or math.isinf(val):
                                    val = 0.0
                            except ValueError:
                                val = 0.0
                            coverm_data[mid_norm][s][metric] = val
        else:
            consolidated_files = sorted([
                os.path.join(cov_dir, f) for f in os.listdir(cov_dir)
                if f.endswith(".tsv") or f.endswith(".txt")
            ])
            if not consolidated_files:
                if not args.dry_run:
                    print(f"[ERROR] No CoverM abundance tables found in {cov_dir}", file=sys.stderr)
                    sys.exit(3)
            for cpath in consolidated_files:
                source_files.append(cpath)
                source_fps[cpath] = compute_sha256(cpath)
                with open(cpath, "r", encoding="utf-8") as f:
                    reader = csv.reader(f, delimiter="\t")
                    header = next(reader, None)
                    if not header or len(header) < 2:
                        continue
                    col_specs = []
                    for col_name in header[1:]:
                        col_s = None
                        for s in samples:
                            if s.lower() in col_name.lower():
                                col_s = s
                                break
                        norm_col = col_name.lower()
                        if col_s:
                            norm_col = norm_col.replace(col_s.lower(), "")
                        norm_col = re.sub(r"[\s\-_]+", "_", norm_col).strip("_").replace("(%)", "").strip("_")
                        if "relative_abundance" in norm_col:
                            canonical = "relative_abundance"
                        elif "covered_fraction" in norm_col:
                            canonical = "covered_fraction"
                        elif "mean" in norm_col and "coverage" in norm_col:
                            canonical = "mean_coverage"
                        elif "trimmed_mean" in norm_col:
                            canonical = "trimmed_mean"
                        elif "covered_bases" in norm_col:
                            canonical = "covered_bases"
                        elif "reads_per_base" in norm_col:
                            canonical = "reads_per_base"
                        elif norm_col in ALLOWED_COVERM_METRICS:
                            canonical = norm_col
                        else:
                            canonical = "relative_abundance"
                        col_specs.append((col_s, canonical))

                    for row in reader:
                        if not row or not row[0].strip():
                            continue
                        mid = row[0].strip()
                        mid_norm = re.sub(r"\.(fa|fasta|fna)$", "", mid, flags=re.IGNORECASE)
                        if mid_norm not in coverm_data:
                            coverm_data[mid_norm] = {}
                        for (col_s, metric), val_str in zip(col_specs, row[1:]):
                            if not col_s:
                                continue
                            if col_s not in coverm_data[mid_norm]:
                                coverm_data[mid_norm][col_s] = {}
                            try:
                                val = float(val_str.strip())
                                if math.isnan(val) or math.isinf(val):
                                    val = 0.0
                            except ValueError:
                                val = 0.0
                            coverm_data[mid_norm][col_s][metric] = val

    if args.dry_run:
        print(f"MAG summary      : would run ({len(base_mags)} base MAGs)")
        sys.exit(0)

    # Build rows
    records = []
    seen_keys = set()

    for mid in base_mags:
        comp, cont = checkm2_data.get(mid, ("NA", "NA"))
        if args.filtering_enabled:
            retained_bool, reason = filter_data.get(mid, (False, "not_evaluated"))
            retained_str = "true" if retained_bool else "false"
            filter_reason = reason
        else:
            retained_bool = True
            retained_str = "NA"
            filter_reason = "NA"

        # Check if excluded by filtering
        if args.filtering_enabled and not retained_bool:
            # Excluded MAG: must not appear in GTDB-Tk or CoverM
            status_tokens = ["quality-reported", "excluded"]
            status = build_mag_status(status_tokens)
            key = (mid, "NA", "NA")
            if key in seen_keys:
                print(f"[ERROR] Duplicate key in MAG summary: {key}", file=sys.stderr)
                sys.exit(1)
            seen_keys.add(key)
            records.append({
                "mag_id": mid,
                "sample_id": "NA",
                "completeness": comp if comp != "NA" else "NA",
                "contamination": cont if cont != "NA" else "NA",
                "retained": retained_str,
                "filter_reason": sanitize_tsv(filter_reason),
                "gtdb_classification": "NA",
                "gtdb_domain": "NA",
                "gtdb_class": "NA",
                "gtdb_order": "NA",
                "gtdb_family": "NA",
                "gtdb_genus": "NA",
                "gtdb_species": "NA",
                "classification_method": "NA",
                "closest_reference": "NA",
                "ani": "NA",
                "af": "NA",
                "coverm_metric": "NA",
                "coverm_value": "NA",
                "status": status,
            })
            continue

        # Retained MAG (or filtering not active)
        base_status_tokens = ["quality-reported"]
        if args.filtering_enabled:
            base_status_tokens.append("retained")
        elif not (args.gtdbtk_enabled or args.coverm_enabled):
            base_status_tokens.append("checkm2-only")

        # GTDB-Tk details
        gtdb_row = gtdb_data.get(mid, None)
        if args.gtdbtk_enabled:
            if gtdb_row and gtdb_row["classification"] not in ("Unclassified", "NA", ""):
                base_status_tokens.append("classified")
                gtdb_clf = gtdb_row["classification"]
                gtdb_lineage = gtdb_row["lineage"]
                clf_method = gtdb_row["method"]
                closest_ref = gtdb_row["closest_reference"]
                ani = gtdb_row["ani"]
                af = gtdb_row["af"]
            else:
                base_status_tokens.append("unclassified")
                gtdb_clf = "Unclassified" if gtdb_row else "NA"
                gtdb_lineage = parse_gtdb_lineage("NA")
                clf_method = "NA"
                closest_ref = "NA"
                ani = "NA"
                af = "NA"
        else:
            gtdb_clf = "NA"
            gtdb_lineage = parse_gtdb_lineage("NA")
            clf_method = "NA"
            closest_ref = "NA"
            ani = "NA"
            af = "NA"

        if args.coverm_enabled:
            configured_methods = [m.strip() for m in args.coverm_methods.split(",") if m.strip()] if args.coverm_methods else ["relative_abundance"]
            target_samples = samples if samples else (sorted(coverm_data[mid].keys()) if coverm_data.get(mid) else [])
            if not target_samples:
                target_samples = ["NA"]

            for s in target_samples:
                for metric in configured_methods:
                    val = coverm_data.get(mid, {}).get(s, {}).get(metric, None) if s != "NA" else None
                    if val is not None:
                        val_str = f"{val:.6g}"
                        metric_status = "coverm-reported" if val > 0 else "coverm-zero"
                    else:
                        val_str = "NA"
                        metric_status = "coverm-missing"
                    status = build_mag_status(base_status_tokens + [metric_status])
                    key = (mid, s, metric)
                    if key in seen_keys:
                        print(f"[ERROR] Duplicate key in MAG summary: {key}", file=sys.stderr)
                        sys.exit(1)
                    seen_keys.add(key)
                    records.append({
                        "mag_id": mid,
                        "sample_id": s,
                        "completeness": comp if comp != "NA" else "NA",
                        "contamination": cont if cont != "NA" else "NA",
                        "retained": retained_str,
                        "filter_reason": sanitize_tsv(filter_reason),
                        "gtdb_classification": sanitize_tsv(gtdb_clf),
                        "gtdb_domain": sanitize_tsv(gtdb_lineage["domain"]),
                        "gtdb_class": sanitize_tsv(gtdb_lineage["class"]),
                        "gtdb_order": sanitize_tsv(gtdb_lineage["order"]),
                        "gtdb_family": sanitize_tsv(gtdb_lineage["family"]),
                        "gtdb_genus": sanitize_tsv(gtdb_lineage["genus"]),
                        "gtdb_species": sanitize_tsv(gtdb_lineage["species"]),
                        "classification_method": sanitize_tsv(clf_method),
                        "closest_reference": sanitize_tsv(closest_ref),
                        "ani": ani,
                        "af": af,
                        "coverm_metric": sanitize_tsv(metric),
                        "coverm_value": val_str,
                        "status": status,
                    })
        else:
            # CoverM disabled or no coverm rows for this MAG
            status = build_mag_status(base_status_tokens)
            key = (mid, "NA", "NA")
            if key in seen_keys:
                print(f"[ERROR] Duplicate key in MAG summary: {key}", file=sys.stderr)
                sys.exit(1)
            seen_keys.add(key)
            records.append({
                "mag_id": mid,
                "sample_id": "NA",
                "completeness": comp if comp != "NA" else "NA",
                "contamination": cont if cont != "NA" else "NA",
                "retained": retained_str,
                "filter_reason": sanitize_tsv(filter_reason),
                "gtdb_classification": sanitize_tsv(gtdb_clf),
                "gtdb_domain": sanitize_tsv(gtdb_lineage["domain"]),
                "gtdb_class": sanitize_tsv(gtdb_lineage["class"]),
                "gtdb_order": sanitize_tsv(gtdb_lineage["order"]),
                "gtdb_family": sanitize_tsv(gtdb_lineage["family"]),
                "gtdb_genus": sanitize_tsv(gtdb_lineage["genus"]),
                "gtdb_species": sanitize_tsv(gtdb_lineage["species"]),
                "classification_method": sanitize_tsv(clf_method),
                "closest_reference": sanitize_tsv(closest_ref),
                "ani": ani,
                "af": af,
                "coverm_metric": "NA",
                "coverm_value": "NA",
                "status": status,
            })

    # Sort deterministically: mag_id, sample_id, coverm_metric
    records.sort(key=lambda r: (r["mag_id"], r["sample_id"], r["coverm_metric"]))

    # Write atomically
    out_dir = args.outdir
    os.makedirs(out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out_dir, prefix=f".build.mag.{os.getpid()}_") as bdir:
        tmp_tsv = os.path.join(bdir, "mag_summary.tsv")
        headers = [
            "mag_id",
            "sample_id",
            "completeness",
            "contamination",
            "retained",
            "filter_reason",
            "gtdb_classification",
            "gtdb_domain",
            "gtdb_class",
            "gtdb_order",
            "gtdb_family",
            "gtdb_genus",
            "gtdb_species",
            "classification_method",
            "closest_reference",
            "ani",
            "af",
            "coverm_metric",
            "coverm_value",
            "status",
        ]
        with open(tmp_tsv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for r in records:
                writer.writerow(r)

        summary_fp = compute_sha256(tmp_tsv)

        # Write metadata JSON
        meta = {
            "schema_version": SCHEMA_VERSION,
            "summary_type": "mag",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "pipeline_version": args.pipeline_version,
            "row_count": len(records),
            "source_files": sorted(list(set(source_files))),
            "source_fingerprints": source_fps,
            "summary_fingerprint": summary_fp,
            "status": "available",
            "bins_manifest": args.bins_dir,
            "checkm2_report": args.checkm2_report if args.checkm2_enabled else "NA",
            "filtering_manifest": args.filter_manifest if args.filtering_enabled else "NA",
            "gtdbtk_outputs": [f for f in source_files if "gtdbtk" in f],
            "coverm_output": [f for f in source_files if "abundance" in f],
            "min_completeness": args.min_completeness,
            "max_contamination": args.max_contamination,
            "coverm_methods": args.coverm_methods if args.coverm_enabled else "NA",
        }
        tmp_meta = os.path.join(bdir, "mag_summary.meta.json")
        with open(tmp_meta, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)

        # Atomic publish
        final_tsv = os.path.join(out_dir, "mag_summary.tsv")
        final_meta = os.path.join(out_dir, "mag_summary.meta.json")
        atomic_publish_file(tmp_tsv, final_tsv)
        atomic_publish_file(tmp_meta, final_meta)

    print(f"[MAG Summary] Published {len(records)} row(s) to {final_tsv}")
    sys.exit(0)


# ------------------------------------------------------------------------------
# Run Overview Subcommand
# ------------------------------------------------------------------------------
KNOWN_STAGES = [
    "01_fastp",
    "03_kraken2",
    "02_megahit",
    "04_mapping",
    "05a_vamb",
    "05b_metabinner",
    "05c_comebin",
    "06_refinement",
    "07_reassembly",
    "08_checkm2",
    "09_filter_mags",
    "10_gtdbtk",
    "11_coverm",
    "12a_taxonomy_summary",
    "12b_mag_summary",
    "12c_run_overview",
]


def run_overview(args):
    # Always applicable if preflight completed
    out_dir = args.outdir
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    seen_keys = set()

    def add_row(section, key, value, source, status):
        ukey = (section, key)
        if ukey in seen_keys:
            return
        seen_keys.add(ukey)
        rows.append({
            "section": section,
            "key": key,
            "value": sanitize_tsv(value),
            "source": source,
            "status": status,
        })

    # 1. Run
    add_row("run", "run_id", args.run_id or "run_default", "configuration", "available")
    add_row("run", "pipeline_version", args.pipeline_version, "configuration", "available")
    add_row("run", "profile", args.profile, "configuration", "available")
    add_row("run", "outdir", args.pipeline_outdir, "configuration", "available")
    add_row("run", "execution_mode", args.execution_mode, "configuration", "available")
    add_row("run", "dry_run", str(args.dry_run).lower(), "configuration", "available")
    add_row("run", "start_time", args.start_time or "unknown", "configuration", "available" if args.start_time else "unknown")
    add_row("run", "end_time", args.end_time or "pending", "configuration", "available" if args.end_time else "pending")
    add_row("run", "final_status", args.final_status or "completed", "run_status", "available")

    # 2. Inputs
    add_row("inputs", "reads_source", args.reads_source or "raw", "configuration", "available")
    add_row("inputs", "sample_manifest", args.sample_manifest or "not-applicable", "sample_manifest", "available" if args.sample_manifest and os.path.isfile(args.sample_manifest) else "not-applicable")
    samples = parse_sample_manifest(args.sample_manifest)
    add_row("inputs", "sample_count", str(len(samples)), "sample_manifest", "available" if samples else "not-applicable")
    add_row("inputs", "sample_manifest_fingerprint", compute_sha256(args.sample_manifest) if args.sample_manifest else "not-applicable", "sample_manifest", "available" if args.sample_manifest and os.path.isfile(args.sample_manifest) else "not-applicable")
    add_row("inputs", "assembly_source", args.assembly_source or "megahit", "configuration", "available")
    add_row("inputs", "assembly_path", args.assembly_path or "not-applicable", "configuration", "available" if args.assembly_path and os.path.isfile(args.assembly_path) else "not-applicable")
    add_row("inputs", "bins_source_requested", args.bins_source_requested or "not-applicable", "configuration", "available" if args.bins_source_requested else "not-applicable")
    add_row("inputs", "bins_source_resolved", args.bins_source_resolved or "not-applicable", "configuration", "available" if args.bins_source_resolved else "not-applicable")
    add_row("inputs", "external_bins_fingerprint", args.external_bins_fingerprint or "not-applicable", "configuration", "available" if args.external_bins_fingerprint else "not-applicable")

    # 3. Configuration
    add_row("configuration", "threads", str(args.threads), "configuration", "available")
    add_row("configuration", "min_completeness", str(args.min_completeness), "configuration", "available")
    add_row("configuration", "max_contamination", str(args.max_contamination), "configuration", "available")
    add_row("configuration", "comebin_device", args.comebin_device or "cpu", "configuration", "available")
    add_row("configuration", "enabled_stages", args.enabled_stages or "all", "configuration", "available")
    add_row("configuration", "refinement_mode", args.refinement_mode or "auto", "configuration", "available")
    add_row("configuration", "reassembly_enabled", str(args.reassembly_enabled).lower(), "configuration", "available")

    # 4. Stages
    plan_actions = {}
    if args.resume_plan and os.path.isfile(args.resume_plan):
        with open(args.resume_plan, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row or not any(row):
                    continue
                if row[0].strip().lower() in ("stage_id", "stage"):
                    continue
                stg = row[0].strip()
                if len(row) >= 4:
                    plan_actions[stg] = row[3].strip()
                elif len(row) >= 2:
                    plan_actions[stg] = row[1].strip()

    status_map = {}
    duration_map = {}
    if getattr(args, "stage_status", None) and os.path.isfile(args.stage_status):
        with open(args.stage_status, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if header:
                stg_idx = header.index("stage_id") if "stage_id" in header else 2
                stat_idx = header.index("final_status") if "final_status" in header else 7
                dur_idx = header.index("wall_time_seconds") if "wall_time_seconds" in header else 10
                run_idx = header.index("run_id") if "run_id" in header else 1
                for row in reader:
                    if not row or not any(row):
                        continue
                    if len(row) > run_idx and getattr(args, "run_id", None) and row[run_idx].strip() != args.run_id:
                        continue
                    if len(row) > stg_idx:
                        stg = row[stg_idx].strip()
                        if len(row) > stat_idx:
                            status_map[stg] = row[stat_idx].strip()
                        if len(row) > dur_idx and row[dur_idx].strip() not in ("NA", ""):
                            duration_map[stg] = row[dur_idx].strip()
    elif args.run_status and os.path.isfile(args.run_status):
        with open(args.run_status, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row or not any(row):
                    continue
                if row[0].strip().lower() in ("stage_id", "stage", "stage_internal_id"):
                    continue
                stg = row[0].strip()
                if len(row) >= 4:
                    status_map[stg] = row[3].strip()
                elif len(row) >= 2:
                    status_map[stg] = row[1].strip()
                if len(row) >= 5:
                    duration_map[stg] = row[4].strip()

    chk_dir = args.checkpoint_dir or os.path.join(args.pipeline_outdir, ".checkpoints")
    for stg in KNOWN_STAGES:
        act = plan_actions.get(stg, "not-applicable")
        st = status_map.get(stg, "not-applicable" if act == "not-applicable" else "pending")
        chk_file = os.path.join(chk_dir, f"{stg}.done")
        chk_exists = os.path.isfile(chk_file)
        chk_status = "available" if chk_exists else ("not-applicable" if act == "not-applicable" else "unavailable")

        add_row("stages", f"{stg}.action", act, "resume_plan", "available" if act != "not-applicable" else "not-applicable")
        add_row("stages", f"{stg}.status", st, "run_status", "available" if st not in ("not-applicable", "failed") else ("failed" if st == "failed" else "not-applicable"))
        if stg in duration_map:
            add_row("stages", f"{stg}.duration_seconds", duration_map[stg], "run_status", "available")
        add_row("stages", f"{stg}.checkpoint", chk_file if chk_exists else "none", "checkpoint", chk_status)
        add_row("stages", f"{stg}.validation_status", "passed" if chk_exists else ("not-applicable" if act == "not-applicable" else "pending"), "checkpoint", chk_status)

    # 5. Outputs
    tax_tsv = os.path.join(out_dir, "taxonomy_summary.tsv")
    mag_tsv = os.path.join(out_dir, "mag_summary.tsv")
    run_tsv = os.path.join(out_dir, "run_overview.tsv")

    tax_stat = args.taxonomy_summary_status or ("available" if os.path.isfile(tax_tsv) else "not-applicable")
    mag_stat = args.mag_summary_status or ("available" if os.path.isfile(mag_tsv) else "not-applicable")

    add_row("outputs", "taxonomy_summary.path", tax_tsv if tax_stat == "available" else "not-applicable", "summary_stage", tax_stat)
    add_row("outputs", "taxonomy_summary.status", tax_stat, "summary_stage", tax_stat)
    add_row("outputs", "mag_summary.path", mag_tsv if mag_stat == "available" else "not-applicable", "summary_stage", mag_stat)
    add_row("outputs", "mag_summary.status", mag_stat, "summary_stage", mag_stat)
    add_row("outputs", "run_overview.path", run_tsv, "summary_stage", "available")
    add_row("outputs", "run_overview.status", "available", "summary_stage", "available")

    # Manifests & Tracking
    run_man = getattr(args, "run_manifest", "")
    stg_stat = getattr(args, "stage_status", "")
    if run_man:
        add_row("manifests", "run_manifest.path", run_man, "configuration", "available" if os.path.isfile(run_man) else "not-applicable")
    if stg_stat:
        add_row("manifests", "stage_status.path", stg_stat, "configuration", "available" if os.path.isfile(stg_stat) else "not-applicable")

    # 6. Software Inventory
    sw_inv = args.software_inventory
    if sw_inv and os.path.isfile(sw_inv) and os.path.getsize(sw_inv) > 0:
        entries = 0
        unk_count = 0
        unav_count = 0
        with open(sw_inv, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            for row in reader:
                if not row or not any(row):
                    continue
                entries += 1
                if len(row) >= 15:
                    if row[14].strip() == "unknown":
                        unk_count += 1
                    elif row[14].strip() == "unavailable":
                        unav_count += 1
        add_row("software", "software_inventory.path", sw_inv, "software_inventory", "available")
        add_row("software", "software_inventory.entry_count", str(entries), "software_inventory", "available")
        add_row("software", "software_inventory.unknown_count", str(unk_count), "software_inventory", "available")
        add_row("software", "software_inventory.unavailable_count", str(unav_count), "software_inventory", "available")
    else:
        add_row("software", "software_inventory.path", sw_inv or "none", "software_inventory", "unavailable")
        add_row("software", "software_inventory.entry_count", "0", "software_inventory", "unavailable")
        add_row("software", "software_inventory.unknown_count", "0", "software_inventory", "unavailable")
        add_row("software", "software_inventory.unavailable_count", "0", "software_inventory", "unavailable")

    # 7. Database Inventory
    db_inv = args.database_inventory
    if db_inv and os.path.isfile(db_inv) and os.path.getsize(db_inv) > 0:
        entries = 0
        unk_count = 0
        unav_count = 0
        with open(db_inv, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            for row in reader:
                if not row or not any(row):
                    continue
                entries += 1
                if len(row) >= 11:
                    if row[10].strip() == "unknown":
                        unk_count += 1
                    elif row[10].strip() == "unavailable":
                        unav_count += 1
        add_row("databases", "database_inventory.path", db_inv, "database_inventory", "available")
        add_row("databases", "database_inventory.entry_count", str(entries), "database_inventory", "available")
        add_row("databases", "database_inventory.unknown_count", str(unk_count), "database_inventory", "available")
        add_row("databases", "database_inventory.unavailable_count", str(unav_count), "database_inventory", "available")
    else:
        add_row("databases", "database_inventory.path", db_inv or "none", "database_inventory", "unavailable")
        add_row("databases", "database_inventory.entry_count", "0", "database_inventory", "unavailable")
        add_row("databases", "database_inventory.unknown_count", "0", "database_inventory", "unavailable")
        add_row("databases", "database_inventory.unavailable_count", "0", "database_inventory", "unavailable")

    # 8. Warnings & Errors
    warn_count = args.warnings_count or 0
    err_count = args.errors_count or 0
    log_ref = args.log_reference or os.path.join(args.pipeline_outdir, "logs")
    add_row("warnings", "warnings.count", str(warn_count), "summary_stage", "available")
    add_row("warnings", "warnings.log_reference", log_ref, "summary_stage", "available")
    add_row("errors", "errors.count", str(err_count), "summary_stage", "available")
    add_row("errors", "errors.log_reference", log_ref, "summary_stage", "available")
    last_err_path = os.path.join(log_ref, "last_error.json")
    if os.path.isfile(last_err_path):
        try:
            with open(last_err_path, "r", encoding="utf-8") as lef:
                last_err = json.load(lef)
            add_row("errors", "errors.error_id", str(last_err.get("error_id", "none")), "last_error", "available")
            add_row("errors", "errors.failed_stage", str(last_err.get("stage_id", "none")), "last_error", "available")
            add_row("errors", "errors.exit_code", str(last_err.get("exit_code", "none")), "last_error", "available")
            add_row("errors", "errors.error_report", last_err_path, "last_error", "available")
        except Exception:
            pass

    if args.dry_run:
        print(f"Run overview     : would run ({len(rows)} entries previewed)")
        sys.exit(0)

    # Write atomically
    with tempfile.TemporaryDirectory(dir=out_dir, prefix=f".build.overview.{os.getpid()}_") as bdir:
        tmp_tsv = os.path.join(bdir, "run_overview.tsv")
        headers = ["section", "key", "value", "source", "status"]
        with open(tmp_tsv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            for r in rows:
                writer.writerow(r)

        summary_fp = compute_sha256(tmp_tsv)

        meta = {
            "schema_version": SCHEMA_VERSION,
            "summary_type": "run-overview",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "pipeline_version": args.pipeline_version,
            "row_count": len(rows),
            "source_files": [
                f for f in (args.resume_plan, args.run_status, getattr(args, "run_manifest", None), getattr(args, "stage_status", None), args.software_inventory, args.database_inventory)
                if f and os.path.isfile(f)
            ],
            "source_fingerprints": {
                f: compute_sha256(f)
                for f in (args.resume_plan, args.run_status, getattr(args, "run_manifest", None), getattr(args, "stage_status", None), args.software_inventory, args.database_inventory)
                if f and os.path.isfile(f)
            },
            "summary_fingerprint": summary_fp,
            "status": "available" if args.final_status != "failed" else "failed",
            "resume_plan_fingerprint": compute_sha256(args.resume_plan),
            "run_status_fingerprint": compute_sha256(args.run_status),
            "software_inventory_fingerprint": compute_sha256(args.software_inventory),
            "database_inventory_fingerprint": compute_sha256(args.database_inventory),
        }
        tmp_meta = os.path.join(bdir, "run_overview.meta.json")
        with open(tmp_meta, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)

        final_tsv = os.path.join(out_dir, "run_overview.tsv")
        final_meta = os.path.join(out_dir, "run_overview.meta.json")
        atomic_publish_file(tmp_tsv, final_tsv)
        atomic_publish_file(tmp_meta, final_meta)

    print(f"[Run Overview] Published {len(rows)} row(s) to {final_tsv}")
    sys.exit(0)


# ------------------------------------------------------------------------------
# Main CLI Parser
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate modular metagenomics pipeline summaries.")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # 1. Taxonomy
    p_tax = subparsers.add_parser("taxonomy", help="Generate taxonomy_summary.tsv")
    p_tax.add_argument("--outdir", required=True, help="Directory OUTDIR/12_summary")
    p_tax.add_argument("--pipeline-version", default="2.0.0")
    p_tax.add_argument("--sample-manifest", required=True)
    p_tax.add_argument("--kraken2-dir", default="")
    p_tax.add_argument("--kraken2-enabled", action="store_true", default=False)
    p_tax.add_argument("--bracken-enabled", action="store_true", default=False)
    p_tax.add_argument("--bracken-level", default="S")
    p_tax.add_argument("--dry-run", action="store_true", default=False)

    # 2. MAG
    p_mag = subparsers.add_parser("mag", help="Generate mag_summary.tsv")
    p_mag.add_argument("--outdir", required=True, help="Directory OUTDIR/12_summary")
    p_mag.add_argument("--pipeline-version", default="2.0.0")
    p_mag.add_argument("--sample-manifest", default="")
    p_mag.add_argument("--bins-dir", default="")
    p_mag.add_argument("--checkm2-report", default="")
    p_mag.add_argument("--checkm2-enabled", action="store_true", default=False)
    p_mag.add_argument("--filter-manifest", default="")
    p_mag.add_argument("--filtering-enabled", action="store_true", default=False)
    p_mag.add_argument("--min-completeness", default="70")
    p_mag.add_argument("--max-contamination", default="10")
    p_mag.add_argument("--gtdbtk-dir", default="")
    p_mag.add_argument("--gtdbtk-enabled", action="store_true", default=False)
    p_mag.add_argument("--coverm-dir", default="")
    p_mag.add_argument("--coverm-enabled", action="store_true", default=False)
    p_mag.add_argument("--coverm-methods", default="relative_abundance")
    p_mag.add_argument("--dry-run", action="store_true", default=False)

    # 3. Overview
    p_ov = subparsers.add_parser("overview", help="Generate run_overview.tsv")
    p_ov.add_argument("--outdir", required=True, help="Directory OUTDIR/12_summary")
    p_ov.add_argument("--pipeline-outdir", required=True)
    p_ov.add_argument("--pipeline-version", default="2.0.0")
    p_ov.add_argument("--run-id", default="")
    p_ov.add_argument("--profile", default="custom")
    p_ov.add_argument("--execution-mode", default="apptainer")
    p_ov.add_argument("--start-time", default="")
    p_ov.add_argument("--end-time", default="")
    p_ov.add_argument("--final-status", default="completed")
    p_ov.add_argument("--resume-plan", default="")
    p_ov.add_argument("--run-status", default="")
    p_ov.add_argument("--run-manifest", default="")
    p_ov.add_argument("--stage-status", default="")
    p_ov.add_argument("--sample-manifest", default="")
    p_ov.add_argument("--reads-source", default="raw")
    p_ov.add_argument("--assembly-path", default="")
    p_ov.add_argument("--assembly-source", default="megahit")
    p_ov.add_argument("--bins-source-requested", default="")
    p_ov.add_argument("--bins-source-resolved", default="")
    p_ov.add_argument("--external-bins-fingerprint", default="")
    p_ov.add_argument("--threads", default="8")
    p_ov.add_argument("--min-completeness", default="70")
    p_ov.add_argument("--max-contamination", default="10")
    p_ov.add_argument("--comebin-device", default="cpu")
    p_ov.add_argument("--enabled-stages", default="")
    p_ov.add_argument("--refinement-mode", default="auto")
    p_ov.add_argument("--reassembly-enabled", default="false")
    p_ov.add_argument("--checkpoint-dir", default="")
    p_ov.add_argument("--taxonomy-summary-status", default="")
    p_ov.add_argument("--mag-summary-status", default="")
    p_ov.add_argument("--software-inventory", default="")
    p_ov.add_argument("--database-inventory", default="")
    p_ov.add_argument("--warnings-count", type=int, default=0)
    p_ov.add_argument("--errors-count", type=int, default=0)
    p_ov.add_argument("--log-reference", default="")
    p_ov.add_argument("--dry-run", action="store_true", default=False)

    args = parser.parse_args()
    if args.subcommand == "taxonomy":
        run_taxonomy(args)
    elif args.subcommand == "mag":
        run_mag(args)
    elif args.subcommand == "overview":
        run_overview(args)
    else:
        sys.exit(2)


if __name__ == "__main__":
    main()
