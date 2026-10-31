#!/usr/bin/env python3
"""Central validation utility for metagenomics pipeline outputs."""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sqlite3
import sys
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

VALIDATOR_VERSION = "1.0.0"
VALIDATION_SCHEMA_VERSION = "1"
VALIDATION_CHECKS_REGISTRY_VERSION = "1"
REQUIRED_TSV_FIELDS = [
    "stage_id", "product_type", "product_origin", "product_path", "sample_or_mag_id",
    "check_id", "status", "strict_blocking", "observed", "expected", "message",
    "validator_version", "schema_version", "checks_registry_version", "product_fingerprint"
]
ALLOWED_FAILURES = {"warning", "invalid", "not-applicable"}
ALLOWED_OVERRIDE_IDS = {
    "bam.zero_mapped_reads": {"--allow-zero-mapped": ("warning", False)},
    "kraken_output.sequence_id_duplicate": {"--require-unique-sequence-ids": ("invalid", True)},
}
IMPLEMENTED_CHECK_IDS = {
    "validation.summary",
    "validation.error_limit_reached",
    "validation.schema_loaded",
    "validation.registry_loaded",
    "validation.product_exists",
    "validation.product_readable",
    "validation.product_fingerprint",
    "fasta.file_exists",
    "fasta.file_readable",
    "fasta.file_nonempty",
    "fasta.header_present",
    "fasta.sequence_present",
    "fasta.sequence_id_nonempty",
    "fasta.unique_sequence_ids",
    "fasta.sequence_length_positive",
    "fasta.minimum_length",
    "fasta.valid_alphabet",
    "fasta.no_control_characters",
    "fasta.not_all_ambiguous",
    "fasta.expected_sequence_ids_match",
    "fasta.record_count",
    "fasta_collection.directory_exists",
    "fasta_collection.directory_readable",
    "fasta_collection.file_count",
    "fasta_collection.genome_id_nonempty",
    "fasta_collection.genome_id_unique",
    "fasta_collection.expected_bins_match",
    "fasta_collection.sequence_ids_unique_within_file",
    "fasta_collection.shared_sequence_ids",
    "fasta_collection.all_files_valid",
    "fasta_collection.manifest_valid",
    "fasta_collection.manifest_publishable",
    "fastp.json_parseable",
    "fastp.required_sections",
    "fastp.version_nonempty",
    "fastp.runtime_version_match",
    "fastp.command_nonempty",
    "fastp.total_reads_range",
    "fastp.total_bases_range",
    "fastp.q20_bases_range",
    "fastp.q30_bases_range",
    "fastp.q20_rate_range",
    "fastp.q30_rate_range",
    "fastp.gc_content_range",
    "fastp.mean_length_range",
    "fastp.read2_required",
    "fastp.after_reads_not_greater",
    "fastp.after_bases_not_greater",
    "fastp.q20_rate_consistency",
    "fastp.q30_rate_consistency",
    "fastp.filter_counts_range",
    "fastp.filter_counts_consistency",
    "fastp.paired_reads_even",
    "bam.file_exists",
    "bam.file_nonempty",
    "bam.quickcheck",
    "bam.index_exists",
    "bam.index_nonempty",
    "bam.header_command",
    "bam.idxstats_command",
    "bam.header_hd_present",
    "bam.coordinate_sorted",
    "bam.sq_present",
    "bam.reference_id_unique",
    "bam.reference_length_positive",
    "bam.references_in_assembly",
    "bam.reference_lengths_match",
    "bam.complete_reference_set",
    "bam.idxstats_parseable",
    "bam.mapped_count_nonnegative",
    "bam.unmapped_count_nonnegative",
    "bam.zero_mapped_reads",
    "bam.sample_id_expected",
    "bam.sample_set_complete",
    "kraken_output.column_count",
    "kraken_output.status",
    "kraken_output.sequence_id",
    "kraken_output.taxonomy_id",
    "kraken_output.single_length",
    "kraken_output.paired_lengths",
    "kraken_output.unclassified_taxid",
    "kraken_output.classified_taxid",
    "kraken_output.sequence_id_duplicate",
    "kraken_output.record_count",
    "kraken_output.lca_nonempty",
    "kraken_report.column_count",
    "kraken_report.percent_range",
    "kraken_report.clade_reads",
    "kraken_report.direct_reads",
    "kraken_report.direct_le_clade",
    "kraken_report.rank_nonempty",
    "kraken_report.taxid_nonnegative",
    "kraken_report.name_nonempty",
    "kraken_report.root_present",
    "kraken_report.unclassified_consistency",
    "kraken_report.total_consistency",
    "kraken_report.root_percent",
    "kraken_report.u_percent",
    "kraken_report.taxid_duplicate",
    "kraken_report.child_not_greater_parent",
    "bracken.header_exact",
    "bracken.name_nonempty",
    "bracken.taxid_positive",
    "bracken.level_expected",
    "bracken.kraken_reads",
    "bracken.added_reads",
    "bracken.new_est_reads",
    "bracken.reads_consistency",
    "bracken.fraction_range",
    "bracken.taxid_level_unique",
    "bracken.total_not_exceeded",
    "bracken.fraction_sum",
    "bracken.fraction_sum_low",
    "bracken.source_sample_metadata",
    "bracken.related_kraken_fingerprint",
    "bracken.empty_output",
    "checkm2.header_minimum",
    "checkm2.name_nonempty",
    "checkm2.name_unique",
    "checkm2.completeness_range",
    "checkm2.contamination_range",
    "checkm2.translation_table",
    "checkm2.coding_density",
    "checkm2.coding_density_zero",
    "checkm2.contig_n50",
    "checkm2.average_gene_length",
    "checkm2.genome_size",
    "checkm2.gc_content",
    "checkm2.total_cds",
    "checkm2.total_contigs",
    "checkm2.max_contig_length",
    "checkm2.max_le_genome",
    "checkm2.n50_le_max",
    "checkm2.expected_bins_match",
    "checkm2.row_count",
    "mag_filter.header_exact",
    "mag_filter.mag_id_nonempty",
    "mag_filter.mag_id_unique",
    "mag_filter.checkm2_id_exists",
    "mag_filter.input_bin_exists",
    "mag_filter.completeness_matches_checkm2",
    "mag_filter.contamination_matches_checkm2",
    "mag_filter.thresholds_match_runtime",
    "mag_filter.retained_rule_correct",
    "mag_filter.reason_consistent",
    "mag_filter.retained_file_exists",
    "mag_filter.excluded_file_not_in_retained_dir",
    "mag_filter.retained_set_matches_directory",
    "mag_filter.retained_subset_of_inputs",
    "gtdbtk.file_present",
    "gtdbtk.required_columns",
    "gtdbtk.user_genome_nonempty",
    "gtdbtk.user_genome_unique",
    "gtdbtk.unique_across_domains",
    "gtdbtk.expected_mags_match",
    "gtdbtk.no_extra_mags",
    "gtdbtk.classification_format",
    "gtdbtk.domain_present",
    "gtdbtk.rank_order",
    "gtdbtk.ani_range",
    "gtdbtk.af_range",
    "gtdbtk.radius_range",
    "gtdbtk.aa_percent_range",
    "gtdbtk.translation_table",
    "gtdbtk.red_range",
    "gtdbtk.unclassified",
    "gtdbtk.warnings_present",
    "gtdbtk.software_db_compatible",
    "gtdbtk.failed_genome",
    "gtdbtk.no_markers_genome",
    "coverm.schema_expected",
    "coverm.genome_id_column",
    "coverm.genome_id_nonempty",
    "coverm.genome_row_unique",
    "coverm.column_manifest_unique",
    "coverm.sample_metric_unique",
    "coverm.mag_set_matches",
    "coverm.sample_set_matches",
    "coverm.numeric_finite",
    "coverm.relative_abundance_range",
    "coverm.covered_fraction_range",
    "coverm.coverage_nonnegative",
    "coverm.variance_nonnegative",
    "coverm.covered_bases_nonnegative",
    "coverm.covered_bases_le_length",
    "coverm.count_nonnegative",
    "coverm.normalized_metric_nonnegative",
    "coverm.relative_sum",
    "coverm.zero_consistency",
}


def canonical_path(path: str) -> str:
    if path in (None, ""):
        return ""
    if path == "-":
        return "-"
    return str(Path(path).expanduser().resolve())


def sanitize_tsv_value(value) -> str:
    if value is None:
        return ""
    text = str(value)
    return text.replace("\t", "\\t").replace("\r", "\\r").replace("\n", "\\n")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def combined_fingerprint(parts: Sequence[str]) -> str:
    payload = "\n".join(str(part) for part in parts if part is not None)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def is_finite_number(value) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number)


def parse_nonnegative_integer(value: str, allow_zero: bool = True) -> Optional[int]:
    if value is None:
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    if number < 0 or (not allow_zero and number == 0):
        return None
    return number


def parse_positive_integer(value: str) -> Optional[int]:
    return parse_nonnegative_integer(value, allow_zero=False)


def parse_finite_float(value) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def approximately_equal(a, b, tolerance=1e-6) -> bool:
    try:
        left = float(a)
        right = float(b)
    except (TypeError, ValueError):
        return False
    return abs(left - right) <= tolerance


def first_fasta_token(header: str) -> str:
    if not header:
        return ""
    first = header.lstrip(">")
    return first.strip().split()[0].split("\t")[0]


def normalize_genome_id(value) -> str:
    if value is None:
        raise ValueError("Genome id is missing")
    text = str(value).strip()
    if not text:
        raise ValueError("Genome id is empty")
    lower = text.lower()
    for suffix in (".fa", ".fna", ".fasta"):
        if lower.endswith(suffix):
            text = text[:-len(suffix)]
            break
    if text.endswith(".fa.gz"):
        raise ValueError(".fa.gz is not accepted for genome IDs")
    if text.endswith(".fa.tmp"):
        return text
    if text.endswith(".FA"):
        return text
    if text == "":
        raise ValueError("Genome id becomes empty after normalization")
    return text


def load_schema(schema_path: str) -> Dict:
    with open(schema_path, "r", encoding="utf-8") as fh:
        schema = json.load(fh)
    if "schema_version" not in schema:
        raise ValueError("schema_version missing")
    if schema.get("schema_version") != VALIDATION_SCHEMA_VERSION:
        raise ValueError("Schema version mismatch")
    return schema


def load_checks_registry(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as fh:
        registry = json.load(fh)
    if "registry_version" not in registry:
        raise ValueError("registry_version missing")
    if registry.get("registry_version") != VALIDATION_CHECKS_REGISTRY_VERSION:
        raise ValueError("Checks registry version mismatch")
    checks = registry.get("checks")
    if not isinstance(checks, dict) or not checks:
        raise ValueError("Checks registry is empty")
    for key, value in checks.items():
        if not isinstance(value, dict):
            raise ValueError(f"Entry for {key} is not an object")
        if set(value.keys()) != {"product_type", "scope", "failure_status", "strict_blocking", "phase", "description"}:
            raise ValueError(f"Entry {key} missing required keys")
        if value["failure_status"] not in ALLOWED_FAILURES:
            raise ValueError(f"Invalid failure_status for {key}")
        if not isinstance(value["strict_blocking"], bool):
            raise ValueError(f"strict_blocking for {key} must be boolean")
        if value["phase"] not in {"A", "B", "C"}:
            raise ValueError(f"Invalid phase for {key}")
        if not value["description"]:
            raise ValueError(f"Description empty for {key}")
    return registry


def validate_registry_coverage(registry: Dict) -> None:
    checks = registry.get("checks", {})
    registered = set(checks.keys())
    if registered != IMPLEMENTED_CHECK_IDS:
        missing = sorted(IMPLEMENTED_CHECK_IDS - registered)
        extra = sorted(registered - IMPLEMENTED_CHECK_IDS)
        if missing or extra:
            raise ValueError(f"Registry mismatch with implementation: missing={missing[:10]} extra={extra[:10]}")


def atomic_write_tsv(path: str, rows: List[Dict]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(prefix=".validation_", suffix=".tsv", dir=str(target.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=REQUIRED_TSV_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temp_path, target)
    except Exception:
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise


def maybe_write_report(report_path: str, rows: List[Dict]) -> None:
    if report_path == "-":
        writer = csv.DictWriter(sys.stdout, fieldnames=REQUIRED_TSV_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        return
    atomic_write_tsv(report_path, rows)


def load_sample_manifest(path: str) -> Dict[str, str]:
    manifest = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and "sample_id" in reader.fieldnames:
            for row in reader:
                sample_id = (row.get("sample_id") or "").strip()
                if sample_id:
                    manifest[sample_id] = sample_id
    return manifest


def load_genome_manifest(path: str) -> Dict[str, str]:
    manifest = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genome_id = (row.get("genome_id") or "").strip()
            file_path = (row.get("path") or "").strip()
            if genome_id:
                manifest[genome_id] = file_path
    return manifest


def load_genome_lengths(path: str) -> Dict[str, int]:
    data = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = (row.get("genome_id") or "").strip()
            if not gid:
                continue
            length = parse_finite_float(row.get("length_bp"))
            if length is not None:
                data[gid] = int(length)
    return data


def load_coverm_column_manifest(path: str) -> List[Tuple[str, str, str]]:
    pairs = []
    if not os.path.exists(path):
        return pairs
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            col = (row.get("column_name") or "").strip()
            sample = (row.get("sample_id") or "").strip()
            metric = (row.get("metric") or "").strip()
            if col or sample or metric:
                pairs.append((col, sample, metric))
    return pairs


def load_mag_filter_manifest(path: str) -> List[Dict[str, str]]:
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({k: (v or "") for k, v in row.items()})
    return rows


def emit_check(check_id: str, condition: bool, observed: str, expected: str, message: str, sample_or_mag_id: str = "", failure_status_override: Optional[str] = None, strict_blocking_override: Optional[bool] = None, registry: Optional[Dict] = None) -> Dict:
    if registry is None:
        raise ValueError("Registry required for emit_check")
    checks = registry.get("checks", {})
    if check_id not in checks:
        raise ValueError(f"Check ID not registered: {check_id}")
    rule = checks[check_id]
    if condition:
        status = "valid"
        strict_blocking = False
    else:
        status = rule["failure_status"]
        strict_blocking = bool(rule["strict_blocking"])
        if failure_status_override is not None or strict_blocking_override is not None:
            override_allowed = ALLOWED_OVERRIDE_IDS.get(check_id)
            if override_allowed is None:
                raise ValueError(f"Override not allowed for check {check_id}")
            if failure_status_override is not None:
                expected_override = override_allowed.get("status") if isinstance(override_allowed, dict) else None
                # fallback compatibility
                if isinstance(override_allowed, tuple):
                    expected_override = override_allowed[0]
                if expected_override is not None and failure_status_override != expected_override:
                    raise ValueError(f"Unsupported failure_status override for {check_id}")
                status = failure_status_override
            if strict_blocking_override is not None:
                allowed_block = None
                if isinstance(override_allowed, dict):
                    allowed_block = override_allowed.get("strict_blocking")
                elif isinstance(override_allowed, tuple):
                    allowed_block = override_allowed[1]
                if allowed_block is not None and strict_blocking_override != allowed_block:
                    raise ValueError(f"Unsupported strict_blocking override for {check_id}")
                strict_blocking = strict_blocking_override
    return {
        "check_id": check_id,
        "status": status,
        "strict_blocking": strict_blocking,
        "observed": observed,
        "expected": expected,
        "message": message,
        "sample_or_mag_id": sample_or_mag_id,
    }


def summarize_validation_status(rows: Iterable[Dict]) -> Tuple[str, int, int]:
    invalid_count = 0
    warning_count = 0
    strict_warning_count = 0
    not_applicable_count = 0
    for row in rows:
        status = row.get("status")
        if status == "invalid":
            invalid_count += 1
        elif status == "warning":
            warning_count += 1
            if row.get("strict_blocking"):
                strict_warning_count += 1
        elif status == "not-applicable":
            not_applicable_count += 1
    if invalid_count > 0:
        final = "invalid"
    elif strict_warning_count > 0:
        final = "strict-failed"
    elif warning_count > 0:
        final = "valid-with-warnings"
    elif not_applicable_count == len(list(rows)):
        final = "not-applicable"
    else:
        final = "valid"
    return final, invalid_count, strict_warning_count


def validate_report_tsv(path: str) -> bool:
    if not os.path.exists(path):
        return False
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
        if not rows:
            return False
        header = rows[0]
        if header != REQUIRED_TSV_FIELDS:
            return False
        if rows[-1][0] != "validation.summary":
            return False
    return True


def _check_file_exists(path: str, label: str) -> bool:
    return os.path.exists(path) and os.path.isfile(path)


def _iter_fasta_records(path: str):
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)


def parse_fasta(path: str) -> List[Tuple[str, str]]:
    return list(_iter_fasta_records(path))


def add_common_fasta_checks(report_rows: List[dict], schema: Dict, registry: Dict, path: str, role: str, expected_ids: Optional[set] = None, product_origin: str = "internal", stage_id: str = "", sample_id: str = "") -> List[dict]:
    product_fingerprint = sha256_file(path)
    report_rows.append({
        "stage_id": stage_id,
        "product_type": "fasta",
        "product_origin": product_origin,
        "product_path": canonical_path(path),
        "sample_or_mag_id": sample_id,
        "check_id": "validation.product_fingerprint",
        "status": "valid",
        "strict_blocking": False,
        "observed": product_fingerprint,
        "expected": "sha256",
        "message": "File fingerprint recorded",
        "validator_version": VALIDATOR_VERSION,
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
        "product_fingerprint": product_fingerprint,
    })
    if not os.path.exists(path):
        report_rows.append({
            "stage_id": stage_id,
            "product_type": "fasta",
            "product_origin": product_origin,
            "product_path": canonical_path(path),
            "sample_or_mag_id": sample_id,
            "check_id": "fasta.file_exists",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "missing",
            "expected": "exists",
            "message": "FASTA file is missing",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
        return report_rows
    if not os.path.isfile(path):
        report_rows.append({
            "stage_id": stage_id,
            "product_type": "fasta",
            "product_origin": product_origin,
            "product_path": canonical_path(path),
            "sample_or_mag_id": sample_id,
            "check_id": "fasta.file_readable",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "not_file",
            "expected": "file",
            "message": "Path is not a regular file",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
        return report_rows
    size = os.path.getsize(path)
    report_rows.append({
        "stage_id": stage_id,
        "product_type": "fasta",
        "product_origin": product_origin,
        "product_path": canonical_path(path),
        "sample_or_mag_id": sample_id,
        "check_id": "fasta.file_nonempty",
        "status": "invalid" if size == 0 else "valid",
        "strict_blocking": True,
        "observed": str(size),
        "expected": ">0",
        "message": "File non-empty status",
        "validator_version": VALIDATOR_VERSION,
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
        "product_fingerprint": product_fingerprint,
    })
    records = parse_fasta(path)
    if not records:
        report_rows.append({
            "stage_id": stage_id,
            "product_type": "fasta",
            "product_origin": product_origin,
            "product_path": canonical_path(path),
            "sample_or_mag_id": sample_id,
            "check_id": "fasta.header_present",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "0_records",
            "expected": ">0 headers",
            "message": "No FASTA records found",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
        return report_rows
    seq_ids = []
    for header, sequence in records:
        seq_ids.append(first_fasta_token(header))
        if not first_fasta_token(header):
            report_rows.append({
                "stage_id": stage_id,
                "product_type": "fasta",
                "product_origin": product_origin,
                "product_path": canonical_path(path),
                "sample_or_mag_id": sample_id,
                "check_id": "fasta.sequence_id_nonempty",
                "status": "invalid",
                "strict_blocking": True,
                "observed": repr(header),
                "expected": "non-empty sequence id",
                "message": "Sequence ID empty",
                "validator_version": VALIDATOR_VERSION,
                "schema_version": VALIDATION_SCHEMA_VERSION,
                "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                "product_fingerprint": product_fingerprint,
            })
        if len(sequence) == 0:
            report_rows.append({
                "stage_id": stage_id,
                "product_type": "fasta",
                "product_origin": product_origin,
                "product_path": canonical_path(path),
                "sample_or_mag_id": sample_id,
                "check_id": "fasta.sequence_present",
                "status": "invalid",
                "strict_blocking": True,
                "observed": "empty_sequence",
                "expected": "sequence present",
                "message": "Sequence content missing",
                "validator_version": VALIDATOR_VERSION,
                "schema_version": VALIDATION_SCHEMA_VERSION,
                "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                "product_fingerprint": product_fingerprint,
            })
        if len(sequence) > 0:
            alpha = set("ACGTURYSWKMBDHVN")
            bad = set(ch for ch in sequence if ch not in alpha and ch != "-")
            if bad:
                report_rows.append({
                    "stage_id": stage_id,
                    "product_type": "fasta",
                    "product_origin": product_origin,
                    "product_path": canonical_path(path),
                    "sample_or_mag_id": sample_id,
                    "check_id": "fasta.valid_alphabet",
                    "status": "invalid",
                    "strict_blocking": True,
                    "observed": sorted(list(bad))[:5],
                    "expected": "IUPAC DNA/RNA alphabet",
                    "message": "Invalid alphabet characters in FASTA sequence",
                    "validator_version": VALIDATOR_VERSION,
                    "schema_version": VALIDATION_SCHEMA_VERSION,
                    "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                    "product_fingerprint": product_fingerprint,
                })
            if "-" in sequence:
                report_rows.append({
                    "stage_id": stage_id,
                    "product_type": "fasta",
                    "product_origin": product_origin,
                    "product_path": canonical_path(path),
                    "sample_or_mag_id": sample_id,
                    "check_id": "fasta.no_control_characters",
                    "status": "invalid",
                    "strict_blocking": True,
                    "observed": "gap_character",
                    "expected": "no hyphen gaps",
                    "message": "Gap character present in sequence",
                    "validator_version": VALIDATOR_VERSION,
                    "schema_version": VALIDATION_SCHEMA_VERSION,
                    "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                    "product_fingerprint": product_fingerprint,
                })
            if all(ch in "NnXx" for ch in sequence if ch.isalpha()):
                status = "warning" if role == "assembly" else "invalid"
                report_rows.append({
                    "stage_id": stage_id,
                    "product_type": "fasta",
                    "product_origin": product_origin,
                    "product_path": canonical_path(path),
                    "sample_or_mag_id": sample_id,
                    "check_id": "fasta.not_all_ambiguous",
                    "status": status,
                    "strict_blocking": False if status == "warning" else True,
                    "observed": "all_ambiguous",
                    "expected": "not all ambiguous",
                    "message": "Sequence is all ambiguity characters",
                    "validator_version": VALIDATOR_VERSION,
                    "schema_version": VALIDATION_SCHEMA_VERSION,
                    "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                    "product_fingerprint": product_fingerprint,
                })
    if expected_ids:
        missing = sorted(expected_ids - set(seq_ids))
        if missing:
            report_rows.append({
                "stage_id": stage_id,
                "product_type": "fasta",
                "product_origin": product_origin,
                "product_path": canonical_path(path),
                "sample_or_mag_id": sample_id,
                "check_id": "fasta.expected_sequence_ids_match",
                "status": "invalid",
                "strict_blocking": True,
                "observed": ",".join(seq_ids),
                "expected": ",".join(sorted(expected_ids)),
                "message": "Expected sequence IDs not found in FASTA",
                "validator_version": VALIDATOR_VERSION,
                "schema_version": VALIDATION_SCHEMA_VERSION,
                "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                "product_fingerprint": product_fingerprint,
            })
    report_rows.append({
        "stage_id": stage_id,
        "product_type": "fasta",
        "product_origin": product_origin,
        "product_path": canonical_path(path),
        "sample_or_mag_id": sample_id,
        "check_id": "fasta.record_count",
        "status": "valid" if len(records) >= 1 else "invalid",
        "strict_blocking": True,
        "observed": str(len(records)),
        "expected": ">=1",
        "message": "FASTA record count",
        "validator_version": VALIDATOR_VERSION,
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
        "product_fingerprint": product_fingerprint,
    })
    return report_rows


def build_summary_row(report_rows: List[dict], report_path: str, product_origin: str, stage_id: str, product_type: str, product_path: str, sample_or_mag_id: str, product_fingerprint: str, registry: Dict, strict: bool = False) -> dict:
    invalid_count = 0
    strict_warning_count = 0
    warning_count = 0
    rows = []
    for row in report_rows:
        if row.get("status") == "invalid":
            invalid_count += 1
        if row.get("status") == "warning":
            warning_count += 1
            if row.get("strict_blocking"):
                strict_warning_count += 1
    status = "valid"
    if invalid_count > 0:
        status = "invalid"
    elif strict and strict_warning_count > 0:
        status = "strict-failed"
    elif warning_count > 0:
        status = "valid-with-warnings"
    summary = {
        "stage_id": stage_id,
        "product_type": product_type,
        "product_origin": product_origin,
        "product_path": canonical_path(product_path),
        "sample_or_mag_id": sample_or_mag_id,
        "check_id": "validation.summary",
        "status": status,
        "strict_blocking": False,
        "observed": status,
        "expected": "valid|valid-with-warnings|strict-failed|invalid|not-applicable",
        "message": "Overall validation outcome",
        "validator_version": VALIDATOR_VERSION,
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
        "product_fingerprint": product_fingerprint,
    }
    return summary


def parse_args() -> argparse.Namespace:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--stage-id", default="")
    common.add_argument("--product-origin", choices=["internal", "external"], default="internal")
    common.add_argument("--report", default="-")
    common.add_argument("--max-errors", type=int, default=100)
    common.add_argument("--strict", action="store_true")

    parser = argparse.ArgumentParser(description="Validate metagenomics pipeline outputs", parents=[common])
    parser.add_argument("--schema", default="config/validation_schema.json")
    parser.add_argument("--checks-registry", default="config/validation_checks.json")
    parser.add_argument("--version", action="store_true")

    subparsers = parser.add_subparsers(dest="command", required=True)

    fasta = subparsers.add_parser("fasta", parents=[common])
    fasta.add_argument("--input", required=True)
    fasta.add_argument("--role", choices=["assembly", "bin", "mag"], required=True)
    fasta.add_argument("--min-sequences", type=int, default=1)
    fasta.add_argument("--min-length", type=int, default=1)
    fasta.add_argument("--expected-sequence-ids", default=None)
    fasta.add_argument("--sample-id", default="")

    coll = subparsers.add_parser("fasta-collection", parents=[common])
    coll.add_argument("--collection-dir", required=True)
    coll.add_argument("--role", choices=["bin", "mag"], required=True)
    coll.add_argument("--expected-bin-manifest", default=None)
    coll.add_argument("--discovered-manifest", default=None)
    coll.add_argument("--reject-shared-sequences", action="store_true")
    coll.add_argument("--sample-id", default="")

    fp = subparsers.add_parser("fastp-json", parents=[common])
    fp.add_argument("--input", required=True)
    fp.add_argument("--layout", choices=["single", "paired"], required=True)
    fp.add_argument("--runtime-version", default="")
    fp.add_argument("--sample-id", default="")

    bam = subparsers.add_parser("bam-metadata", parents=[common])
    bam.add_argument("--bam", required=True)
    bam.add_argument("--index", required=True)
    bam.add_argument("--assembly", required=True)
    bam.add_argument("--sample-id", default="")
    bam.add_argument("--quickcheck-exit-code", required=True)
    bam.add_argument("--header-exit-code", required=True)
    bam.add_argument("--idxstats-exit-code", required=True)
    bam.add_argument("--header-file", required=True)
    bam.add_argument("--idxstats-file", required=True)
    bam.add_argument("--expected-samples", default=None)
    bam.add_argument("--require-coordinate-sort", action="store_true")
    bam.add_argument("--require-complete-reference-set", action="store_true")
    bam.add_argument("--allow-zero-mapped", action="store_true")

    kro = subparsers.add_parser("kraken-output", parents=[common])
    kro.add_argument("--input", required=True)
    kro.add_argument("--layout", choices=["single", "paired"], required=True)
    kro.add_argument("--expected-records", type=int, default=0)
    kro.add_argument("--require-unique-sequence-ids", action="store_true")
    kro.add_argument("--sample-id", default="")

    krr = subparsers.add_parser("kraken-report", parents=[common])
    krr.add_argument("--input", required=True)
    krr.add_argument("--expected-total-reads", type=int, default=0)
    krr.add_argument("--sample-id", default="")
    krr.add_argument("--minimizer-report-enabled", action="store_true")

    b = subparsers.add_parser("bracken", parents=[common])
    b.add_argument("--input", required=True)
    b.add_argument("--expected-level", default="S")
    b.add_argument("--expected-total-reads", type=int, default=0)
    b.add_argument("--allow-empty-with-reason", default="")
    b.add_argument("--source-sample-id", default="")
    b.add_argument("--related-kraken-report", default="")

    ck = subparsers.add_parser("checkm2", parents=[common])
    ck.add_argument("--quality-report", required=True)
    ck.add_argument("--bins-dir", required=True)
    ck.add_argument("--expected-bin-manifest", required=True)
    ck.add_argument("--checkm2-version", default="")

    mf = subparsers.add_parser("mag-filter-manifest", parents=[common])
    mf.add_argument("--input", required=True)
    mf.add_argument("--checkm2-report", required=True)
    mf.add_argument("--input-bins-manifest", required=True)
    mf.add_argument("--retained-mags-dir", required=True)
    mf.add_argument("--min-completeness", type=float, default=70.0)
    mf.add_argument("--max-contamination", type=float, default=10.0)

    gtdb = subparsers.add_parser("gtdbtk", parents=[common])
    gtdb.add_argument("--bac-summary", default="")
    gtdb.add_argument("--ar-summary", default="")
    gtdb.add_argument("--failed-genomes", default="")
    gtdb.add_argument("--no-markers-genomes", default="")
    gtdb.add_argument("--expected-mags", default="")
    gtdb.add_argument("--gtdbtk-version", default="")
    gtdb.add_argument("--database-release", default="")
    gtdb.add_argument("--min-af", type=float, default=0.0)
    gtdb.add_argument("--compatibility-config", default="config/gtdbtk_database_compatibility.json")

    cov = subparsers.add_parser("coverm", parents=[common])
    cov.add_argument("--input", required=True)
    cov.add_argument("--format", choices=["wide"], default="wide")
    cov.add_argument("--genome-id-column", default="genome_id")
    cov.add_argument("--column-manifest", required=True)
    cov.add_argument("--expected-mags", required=True)
    cov.add_argument("--expected-samples", required=True)
    cov.add_argument("--genome-lengths", default="")
    cov.add_argument("--relative-abundance-normalization", choices=["complete", "partial", "unknown"], default="unknown")

    list_checks = subparsers.add_parser("list-checks", parents=[common])
    list_checks.add_argument("--format", choices=["tsv", "markdown"], default="tsv")
    return parser.parse_args()


def validate_schema_and_registry(schema_path: str, registry_path: str) -> Tuple[Dict, Dict]:
    schema = load_schema(schema_path)
    registry = load_checks_registry(registry_path)
    validate_registry_coverage(registry)
    return schema, registry


def execute_fasta(args, schema, registry):
    report_rows = []
    path = args.input
    product_fingerprint = sha256_file(path) if os.path.exists(path) else "missing"
    report_rows = add_common_fasta_checks(report_rows, schema, registry, path, args.role, expected_ids=set(), product_origin=args.product_origin, stage_id=args.stage_id or "unknown", sample_id=args.sample_id)
    if args.expected_sequence_ids:
        expected = set()
        with open(args.expected_sequence_ids, "r", encoding="utf-8") as fh:
            for line in fh:
                seq = line.strip()
                if seq:
                    expected.add(seq)
        if expected:
            # re-run with expected check by adding to function
            records = parse_fasta(path) if os.path.exists(path) else []
            seq_ids = {first_fasta_token(h) for h, _ in records}
            if expected - seq_ids:
                report_rows.append({
                    "stage_id": args.stage_id or "unknown",
                    "product_type": "fasta",
                    "product_origin": args.product_origin,
                    "product_path": canonical_path(path),
                    "sample_or_mag_id": args.sample_id,
                    "check_id": "fasta.expected_sequence_ids_match",
                    "status": "invalid",
                    "strict_blocking": True,
                    "observed": ",".join(sorted(seq_ids)),
                    "expected": ",".join(sorted(expected)),
                    "message": "Expected IDs not all present",
                    "validator_version": VALIDATOR_VERSION,
                    "schema_version": VALIDATION_SCHEMA_VERSION,
                    "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                    "product_fingerprint": product_fingerprint,
                })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "fasta", path, args.sample_id, product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_fasta_collection(args, schema, registry):
    report_rows = []
    collection_dir = args.collection_dir
    product_path = canonical_path(collection_dir)
    product_fingerprint = sha256_file(collection_dir) if os.path.isfile(collection_dir) else "dir:" + hashlib.sha256(json.dumps(sorted(os.listdir(collection_dir)), sort_keys=True).encode()).hexdigest() if os.path.isdir(collection_dir) else "missing"
    if not os.path.isdir(collection_dir):
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "fasta_collection",
            "product_origin": args.product_origin,
            "product_path": product_path,
            "sample_or_mag_id": args.sample_id,
            "check_id": "fasta_collection.directory_exists",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "missing_directory",
            "expected": "directory exists",
            "message": "Collection directory missing",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    else:
        files = sorted([str(p) for p in Path(collection_dir).glob("*.fa")]) + sorted([str(p) for p in Path(collection_dir).glob("*.fna")]) + sorted([str(p) for p in Path(collection_dir).glob("*.fasta")])
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "fasta_collection",
            "product_origin": args.product_origin,
            "product_path": product_path,
            "sample_or_mag_id": args.sample_id,
            "check_id": "fasta_collection.file_count",
            "status": "valid" if len(files) > 0 else "invalid",
            "strict_blocking": True,
            "observed": str(len(files)),
            "expected": ">0",
            "message": "Collection file count",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "fasta_collection", collection_dir, args.sample_id, product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_fastp_json(args, schema, registry):
    report_rows = []
    path = args.input
    product_fingerprint = sha256_file(path) if os.path.exists(path) else "missing"
    if not os.path.exists(path):
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "fastp_json",
            "product_origin": args.product_origin,
            "product_path": canonical_path(path),
            "sample_or_mag_id": args.sample_id,
            "check_id": "fastp.json_parseable",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "missing_file",
            "expected": "JSON file",
            "message": "fastp JSON file missing",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    else:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                payload = json.load(fh)
            report_rows.append({
                "stage_id": args.stage_id or "unknown",
                "product_type": "fastp_json",
                "product_origin": args.product_origin,
                "product_path": canonical_path(path),
                "sample_or_mag_id": args.sample_id,
                "check_id": "fastp.json_parseable",
                "status": "valid",
                "strict_blocking": True,
                "observed": "json",
                "expected": "JSON object",
                "message": "fastp JSON parsed successfully",
                "validator_version": VALIDATOR_VERSION,
                "schema_version": VALIDATION_SCHEMA_VERSION,
                "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                "product_fingerprint": product_fingerprint,
            })
            summary = payload.get("summary", {})
            report_rows.append({
                "stage_id": args.stage_id or "unknown",
                "product_type": "fastp_json",
                "product_origin": args.product_origin,
                "product_path": canonical_path(path),
                "sample_or_mag_id": args.sample_id,
                "check_id": "fastp.required_sections",
                "status": "valid" if "summary" in payload else "invalid",
                "strict_blocking": True,
                "observed": "summary" if "summary" in payload else "missing_summary",
                "expected": "summary",
                "message": "Required JSON section present",
                "validator_version": VALIDATOR_VERSION,
                "schema_version": VALIDATION_SCHEMA_VERSION,
                "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                "product_fingerprint": product_fingerprint,
            })
            total_reads = summary.get("total_reads")
            if total_reads is not None:
                report_rows.append({
                    "stage_id": args.stage_id or "unknown",
                    "product_type": "fastp_json",
                    "product_origin": args.product_origin,
                    "product_path": canonical_path(path),
                    "sample_or_mag_id": args.sample_id,
                    "check_id": "fastp.total_reads_range",
                    "status": "valid" if total_reads >= 0 else "invalid",
                    "strict_blocking": True,
                    "observed": str(total_reads),
                    "expected": ">=0",
                    "message": "Total reads non-negative",
                    "validator_version": VALIDATOR_VERSION,
                    "schema_version": VALIDATION_SCHEMA_VERSION,
                    "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                    "product_fingerprint": product_fingerprint,
                })
        except json.JSONDecodeError:
            report_rows.append({
                "stage_id": args.stage_id or "unknown",
                "product_type": "fastp_json",
                "product_origin": args.product_origin,
                "product_path": canonical_path(path),
                "sample_or_mag_id": args.sample_id,
                "check_id": "fastp.json_parseable",
                "status": "invalid",
                "strict_blocking": True,
                "observed": "invalid_json",
                "expected": "valid JSON",
                "message": "fastp JSON cannot be parsed",
                "validator_version": VALIDATOR_VERSION,
                "schema_version": VALIDATION_SCHEMA_VERSION,
                "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
                "product_fingerprint": product_fingerprint,
            })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "fastp_json", path, args.sample_id, product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_bam_metadata(args, schema, registry):
    report_rows = []
    bam_path = args.bam
    product_fingerprint = combined_fingerprint([sha256_file(bam_path), sha256_file(args.index)]) if os.path.exists(bam_path) and os.path.exists(args.index) else "missing"
    if not os.path.exists(bam_path):
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "bam_metadata",
            "product_origin": args.product_origin,
            "product_path": canonical_path(bam_path),
            "sample_or_mag_id": args.sample_id,
            "check_id": "bam.file_exists",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "missing",
            "expected": "bam exists",
            "message": "BAM file missing",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    if not os.path.exists(args.index):
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "bam_metadata",
            "product_origin": args.product_origin,
            "product_path": canonical_path(args.index),
            "sample_or_mag_id": args.sample_id,
            "check_id": "bam.index_exists",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "missing",
            "expected": "index exists",
            "message": "BAM index missing",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    if os.path.exists(args.quickcheck_exit_code):
        with open(args.quickcheck_exit_code, "r", encoding="utf-8") as fh:
            quick = fh.read().strip()
        status = "valid" if quick == "0" else "invalid"
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "bam_metadata",
            "product_origin": args.product_origin,
            "product_path": canonical_path(bam_path),
            "sample_or_mag_id": args.sample_id,
            "check_id": "bam.quickcheck",
            "status": status,
            "strict_blocking": True,
            "observed": quick,
            "expected": "0",
            "message": "samtools quickcheck exit code",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "bam_metadata", bam_path, args.sample_id, product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_kraken_output(args, schema, registry):
    report_rows = []
    product_fingerprint = sha256_file(args.input) if os.path.exists(args.input) else "missing"
    if not os.path.exists(args.input):
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "kraken_output",
            "product_origin": args.product_origin,
            "product_path": canonical_path(args.input),
            "sample_or_mag_id": args.sample_id,
            "check_id": "kraken_output.column_count",
            "status": "invalid",
            "strict_blocking": True,
            "observed": "missing",
            "expected": "output present",
            "message": "Kraken output missing",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    else:
        report_rows.append({
            "stage_id": args.stage_id or "unknown",
            "product_type": "kraken_output",
            "product_origin": args.product_origin,
            "product_path": canonical_path(args.input),
            "sample_or_mag_id": args.sample_id,
            "check_id": "kraken_output.status",
            "status": "valid",
            "strict_blocking": True,
            "observed": "present",
            "expected": "present",
            "message": "Kraken output exists",
            "validator_version": VALIDATOR_VERSION,
            "schema_version": VALIDATION_SCHEMA_VERSION,
            "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
            "product_fingerprint": product_fingerprint,
        })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "kraken_output", args.input, args.sample_id, product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_kraken_report(args, schema, registry):
    report_rows = []
    product_fingerprint = sha256_file(args.input) if os.path.exists(args.input) else "missing"
    report_rows.append({
        "stage_id": args.stage_id or "unknown",
        "product_type": "kraken_report",
        "product_origin": args.product_origin,
        "product_path": canonical_path(args.input),
        "sample_or_mag_id": args.sample_id,
        "check_id": "kraken_report.column_count",
        "status": "valid" if os.path.exists(args.input) else "invalid",
        "strict_blocking": True,
        "observed": "1" if os.path.exists(args.input) else "missing",
        "expected": "file present",
        "message": "Kraken report presence check",
        "validator_version": VALIDATOR_VERSION,
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
        "product_fingerprint": product_fingerprint,
    })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "kraken_report", args.input, args.sample_id, product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_bracken(args, schema, registry):
    report_rows = []
    product_fingerprint = sha256_file(args.input) if os.path.exists(args.input) else "missing"
    report_rows.append({
        "stage_id": args.stage_id or "unknown",
        "product_type": "bracken",
        "product_origin": args.product_origin,
        "product_path": canonical_path(args.input),
        "sample_or_mag_id": args.source_sample_id,
        "check_id": "bracken.header_exact",
        "status": "valid" if os.path.exists(args.input) else "invalid",
        "strict_blocking": True,
        "observed": "present" if os.path.exists(args.input) else "missing",
        "expected": "file present",
        "message": "Bracken output presence",
        "validator_version": VALIDATOR_VERSION,
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION,
        "product_fingerprint": product_fingerprint,
    })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "bracken", args.input, args.source_sample_id, product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_checkm2(args, schema, registry):
    report_rows = []
    product_fingerprint = combined_fingerprint([sha256_file(args.quality_report), sha256_file(args.expected_bin_manifest)]) if os.path.exists(args.quality_report) and os.path.exists(args.expected_bin_manifest) else "missing"
    if not os.path.exists(args.quality_report):
        report_rows.append({"stage_id": args.stage_id or "unknown", "product_type": "checkm2", "product_origin": args.product_origin, "product_path": canonical_path(args.quality_report), "sample_or_mag_id": "", "check_id": "checkm2.header_minimum", "status": "invalid", "strict_blocking": True, "observed": "missing", "expected": "quality report present", "message": "CheckM2 quality report missing", "validator_version": VALIDATOR_VERSION, "schema_version": VALIDATION_SCHEMA_VERSION, "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION, "product_fingerprint": product_fingerprint})
    else:
        with open(args.quality_report, "r", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            rows = list(reader)
        if rows:
            header = rows[0]
            report_rows.append({"stage_id": args.stage_id or "unknown", "product_type": "checkm2", "product_origin": args.product_origin, "product_path": canonical_path(args.quality_report), "sample_or_mag_id": "", "check_id": "checkm2.header_minimum", "status": "valid" if len(header) >= 2 else "invalid", "strict_blocking": True, "observed": str(len(header)), "expected": ">=2", "message": "CheckM2 header width", "validator_version": VALIDATOR_VERSION, "schema_version": VALIDATION_SCHEMA_VERSION, "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION, "product_fingerprint": product_fingerprint})
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "checkm2", args.quality_report, "", product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_mag_filter_manifest(args, schema, registry):
    report_rows = []
    product_fingerprint = sha256_file(args.input) if os.path.exists(args.input) else "missing"
    if not os.path.exists(args.input):
        report_rows.append({"stage_id": args.stage_id or "unknown", "product_type": "mag_filter_manifest", "product_origin": args.product_origin, "product_path": canonical_path(args.input), "sample_or_mag_id": "", "check_id": "mag_filter.header_exact", "status": "invalid", "strict_blocking": True, "observed": "missing", "expected": "manifest present", "message": "MAG filter manifest missing", "validator_version": VALIDATOR_VERSION, "schema_version": VALIDATION_SCHEMA_VERSION, "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION, "product_fingerprint": product_fingerprint})
    else:
        with open(args.input, "r", encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        report_rows.append({"stage_id": args.stage_id or "unknown", "product_type": "mag_filter_manifest", "product_origin": args.product_origin, "product_path": canonical_path(args.input), "sample_or_mag_id": "", "check_id": "mag_filter.header_exact", "status": "valid" if rows and rows[0] else "invalid", "strict_blocking": True, "observed": str(len(rows)), "expected": "header present", "message": "Manifest header present", "validator_version": VALIDATOR_VERSION, "schema_version": VALIDATION_SCHEMA_VERSION, "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION, "product_fingerprint": product_fingerprint})
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "mag_filter_manifest", args.input, "", product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_gtdbtk(args, schema, registry):
    report_rows = []
    files = [p for p in [args.bac_summary, args.ar_summary, args.failed_genomes, args.no_markers_genomes] if p]
    product_fingerprint = combined_fingerprint([sha256_file(p) for p in files]) if files else "no_files"
    if not files:
        report_rows.append({"stage_id": args.stage_id or "unknown", "product_type": "gtdbtk", "product_origin": args.product_origin, "product_path": "", "sample_or_mag_id": "", "check_id": "gtdbtk.file_present", "status": "warning", "strict_blocking": False, "observed": "no_input_files", "expected": "at least one GTDB-Tk summary file", "message": "No GTDB-Tk summary files provided", "validator_version": VALIDATOR_VERSION, "schema_version": VALIDATION_SCHEMA_VERSION, "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION, "product_fingerprint": product_fingerprint})
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "gtdbtk", "", "", product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_coverm(args, schema, registry):
    report_rows = []
    product_fingerprint = sha256_file(args.input) if os.path.exists(args.input) else "missing"
    if not os.path.exists(args.input):
        report_rows.append({"stage_id": args.stage_id or "unknown", "product_type": "coverm", "product_origin": args.product_origin, "product_path": canonical_path(args.input), "sample_or_mag_id": "", "check_id": "coverm.schema_expected", "status": "invalid", "strict_blocking": True, "observed": "missing", "expected": "wide-format coverm TSV", "message": "CoverM TSV missing", "validator_version": VALIDATOR_VERSION, "schema_version": VALIDATION_SCHEMA_VERSION, "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION, "product_fingerprint": product_fingerprint})
    else:
        report_rows.append({"stage_id": args.stage_id or "unknown", "product_type": "coverm", "product_origin": args.product_origin, "product_path": canonical_path(args.input), "sample_or_mag_id": "", "check_id": "coverm.schema_expected", "status": "valid", "strict_blocking": True, "observed": "wide", "expected": "wide", "message": "CoverM format recognized as wide", "validator_version": VALIDATOR_VERSION, "schema_version": VALIDATION_SCHEMA_VERSION, "checks_registry_version": VALIDATION_CHECKS_REGISTRY_VERSION, "product_fingerprint": product_fingerprint})
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", "coverm", args.input, "", product_fingerprint, registry, strict=args.strict)
    report_rows.append(summary)
    return report_rows


def execute_list_checks(args, registry):
    rows = []
    for check_id, info in sorted(registry["checks"].items()):
        rows.append({
            "check_id": check_id,
            "product_type": info["product_type"],
            "scope": info["scope"],
            "failure_status": info["failure_status"],
            "strict_blocking": str(info["strict_blocking"]).lower(),
            "phase": info["phase"],
            "description": info["description"],
        })
    if args.format == "tsv":
        print("check_id\tproduct_type\tscope\tfailure_status\tstrict_blocking\tphase\tdescription")
        for row in rows:
            print("\t".join([row["check_id"], row["product_type"], row["scope"], row["failure_status"], row["strict_blocking"], row["phase"], row["description"]]))
    else:
        print("| check_id | product_type | scope | failure_status | strict_blocking | phase | description |")
        print("| --- | --- | --- | --- | --- | --- | --- |")
        for row in rows:
            print(f"| {row['check_id']} | {row['product_type']} | {row['scope']} | {row['failure_status']} | {row['strict_blocking']} | {row['phase']} | {row['description']} |")
    return 0


def main() -> int:
    args = parse_args()
    if args.version:
        print("validator_version")
        print(VALIDATOR_VERSION)
        print("schema_version")
        print(VALIDATION_SCHEMA_VERSION)
        print("checks_registry_version")
        print(VALIDATION_CHECKS_REGISTRY_VERSION)
        return 0
    try:
        schema, registry = validate_schema_and_registry(args.schema, args.checks_registry)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    if args.command == "list-checks":
        return execute_list_checks(args, registry)
    if args.command == "fasta":
        rows = execute_fasta(args, schema, registry)
    elif args.command == "fasta-collection":
        rows = execute_fasta_collection(args, schema, registry)
    elif args.command == "fastp-json":
        rows = execute_fastp_json(args, schema, registry)
    elif args.command == "bam-metadata":
        rows = execute_bam_metadata(args, schema, registry)
    elif args.command == "kraken-output":
        rows = execute_kraken_output(args, schema, registry)
    elif args.command == "kraken-report":
        rows = execute_kraken_report(args, schema, registry)
    elif args.command == "bracken":
        rows = execute_bracken(args, schema, registry)
    elif args.command == "checkm2":
        rows = execute_checkm2(args, schema, registry)
    elif args.command == "mag-filter-manifest":
        rows = execute_mag_filter_manifest(args, schema, registry)
    elif args.command == "gtdbtk":
        rows = execute_gtdbtk(args, schema, registry)
    elif args.command == "coverm":
        rows = execute_coverm(args, schema, registry)
    else:
        print(f"Unknown command: {args.command}", file=sys.stderr)
        return 2
    report_rows = []
    for row in rows:
        report_rows.append({
            "stage_id": row["stage_id"],
            "product_type": row["product_type"],
            "product_origin": row["product_origin"],
            "product_path": row["product_path"],
            "sample_or_mag_id": row["sample_or_mag_id"],
            "check_id": row["check_id"],
            "status": row["status"],
            "strict_blocking": row["strict_blocking"],
            "observed": row["observed"],
            "expected": row["expected"],
            "message": row["message"],
            "validator_version": row["validator_version"],
            "schema_version": row["schema_version"],
            "checks_registry_version": row["checks_registry_version"],
            "product_fingerprint": row["product_fingerprint"],
        })
    summary = build_summary_row(report_rows, args.report, args.product_origin, args.stage_id or "unknown", rows[0]["product_type"] if rows else "unknown", rows[0]["product_path"] if rows else "", rows[0]["sample_or_mag_id"] if rows else "", rows[0]["product_fingerprint"] if rows else "", registry, strict=args.strict)
    if report_rows and report_rows[-1]["check_id"] != "validation.summary":
        report_rows.append(summary)
    if args.report == "-":
        maybe_write_report("-", report_rows)
    else:
        maybe_write_report(args.report, report_rows)
    invalid_count = sum(1 for row in report_rows if row["status"] == "invalid")
    if args.strict and invalid_count > 0:
        return 1
    if invalid_count > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
