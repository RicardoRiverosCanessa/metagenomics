#!/usr/bin/env python3
"""Unit and integration tests for MetaNexus execution plan model, renderer, and CLI."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.execution_plan import (
    ALL_STAGE_FIELDS,
    CANONICAL_STAGE_IDS,
    EXIT_BLOCKED_OR_STRICT,
    EXIT_INTERNAL_ERROR,
    EXIT_OK,
    EXIT_USAGE_ERROR,
    PIPELINE_NAME,
    PIPELINE_RELEASE_STATUS,
    PIPELINE_VERSION,
    SCHEMA_VERSION,
    VALID_APPLICABILITY,
    VALID_CHECKPOINT_STATUS,
    VALID_DATABASE_STATUS,
    VALID_INPUT_STATUS,
    VALID_OUTPUT_STATUS,
    VALID_PLANNED_ACTIONS,
    VALID_REQUESTED_ACTIONS,
    VALID_RUNTIME_STATUS,
    VALID_VALIDATION_STATUS,
    ExecutionPlanInternalError,
    ExecutionPlanValidationError,
    assemble_plan_from_tsv,
    canonicalize_plan_data,
    compute_plan_summary,
    is_plan_executable,
    render_json,
    render_text,
    validate_execution_plan,
)

SCRIPT_PATH = REPO_ROOT / "run_metagenomics_pipeline.sh"


def make_valid_stage_record(stage_id: str = "fastp", **kwargs: Any) -> Dict[str, Any]:
    """Helper to create a valid minimal stage record."""
    rec: Dict[str, Any] = {
        "stage_id": stage_id,
        "stage_name": stage_id.upper(),
        "configured_enabled": True,
        "required": True,
        "requested_action": "normal",
        "applicability": "applicable",
        "planned_action": "run",
        "checkpoint_path": f"/tmp/checkpoints/{stage_id}.done",
        "checkpoint_status": "missing",
        "checkpoint_reason": None,
        "input_status": "available",
        "runtime_status": "available",
        "database_status": "not-applicable",
        "output_status": "missing",
        "validation_status": "not-applicable",
        "invalidated_by": [],
        "blocking_stage": None,
        "reason_code": "checkpoint-missing",
        "reason_detail": None,
        "diagnostic_codes": [],
        "dependencies": [],
        "downstream": [],
        "primary_output": f"/tmp/out/{stage_id}/out.txt",
        "primary_log": f"/tmp/out/logs/{stage_id}.log",
    }
    rec.update(kwargs)
    return rec


def make_valid_plan_doc(stages: list = None, profile: str = "full") -> Dict[str, Any]:
    """Helper to create a valid minimal plan document."""
    if stages is None:
        stages = [make_valid_stage_record("fastp")]
    return {
        "schema_version": SCHEMA_VERSION,
        "pipeline": {
            "name": PIPELINE_NAME,
            "version": PIPELINE_VERSION,
            "release_status": PIPELINE_RELEASE_STATUS,
            "profile": profile,
        },
        "plan": stages,
    }


class TestModelSchemaValidation(unittest.TestCase):
    """Test schema validation, field types, enums, and data integrity."""

    def test_minimal_valid_json(self):
        doc = make_valid_plan_doc()
        # Should not raise
        validate_execution_plan(doc)

    def test_schema_version_mismatch(self):
        doc = make_valid_plan_doc()
        doc["schema_version"] = "2"
        with self.assertRaises(ExecutionPlanInternalError):
            validate_execution_plan(doc)

    def test_pipeline_identity_mismatch(self):
        doc = make_valid_plan_doc()
        doc["pipeline"]["name"] = "OtherPipeline"
        with self.assertRaises(ExecutionPlanInternalError):
            validate_execution_plan(doc)

        doc = make_valid_plan_doc()
        doc["pipeline"]["version"] = "9.9.9"
        with self.assertRaises(ExecutionPlanInternalError):
            validate_execution_plan(doc)

        doc = make_valid_plan_doc()
        doc["pipeline"]["release_status"] = "production"
        with self.assertRaises(ExecutionPlanInternalError):
            validate_execution_plan(doc)

    def test_missing_required_top_level_fields(self):
        for field in ["schema_version", "pipeline", "plan"]:
            doc = make_valid_plan_doc()
            del doc[field]
            with self.assertRaises(ExecutionPlanValidationError):
                validate_execution_plan(doc)

    def test_canonical_stage_record_has_24_fields(self):
        self.assertEqual(len(ALL_STAGE_FIELDS), 24)
        for field in ALL_STAGE_FIELDS:
            doc = make_valid_plan_doc([make_valid_stage_record()])
            del doc["plan"][0][field]
            with self.assertRaises(ExecutionPlanValidationError):
                validate_execution_plan(doc)

    def test_closed_enums_respected(self):
        # Planned actions
        for act in VALID_PLANNED_ACTIONS:
            req = "force" if act == "force" else ("restart" if act == "restart" else "normal")
            ckpt = "valid" if act == "reuse" else "missing"
            conf = False if act == "disabled" else True
            appl = "not-applicable" if act == "not-applicable" else "applicable"
            req_flag = False if act in {"disabled", "not-applicable"} else True
            doc = make_valid_plan_doc([
                make_valid_stage_record(
                    planned_action=act,
                    requested_action=req,
                    checkpoint_status=ckpt,
                    configured_enabled=conf,
                    applicability=appl,
                    required=req_flag,
                )
            ])
            validate_execution_plan(doc)
        with self.assertRaises(ExecutionPlanValidationError):
            validate_execution_plan(make_valid_plan_doc([make_valid_stage_record(planned_action="skip")]))

        # Requested actions
        for act in VALID_REQUESTED_ACTIONS:
            doc = make_valid_plan_doc([make_valid_stage_record(requested_action=act)])
            validate_execution_plan(doc)
        with self.assertRaises(ExecutionPlanValidationError):
            validate_execution_plan(make_valid_plan_doc([make_valid_stage_record(requested_action="rerun")]))

        # Applicability
        for appl in VALID_APPLICABILITY:
            req_flag = False if appl == "not-applicable" else True
            doc = make_valid_plan_doc([make_valid_stage_record(applicability=appl, required=req_flag)])
            validate_execution_plan(doc)
        with self.assertRaises(ExecutionPlanValidationError):
            validate_execution_plan(make_valid_plan_doc([make_valid_stage_record(applicability="unknown")]))

        # Checkpoint status
        for st in VALID_CHECKPOINT_STATUS:
            doc = make_valid_plan_doc([make_valid_stage_record(checkpoint_status=st)])
            validate_execution_plan(doc)
        with self.assertRaises(ExecutionPlanValidationError):
            validate_execution_plan(make_valid_plan_doc([make_valid_stage_record(checkpoint_status="bad")]))

        # Input status
        for st in VALID_INPUT_STATUS:
            doc = make_valid_plan_doc([make_valid_stage_record(input_status=st)])
            validate_execution_plan(doc)

        # Runtime status
        for st in VALID_RUNTIME_STATUS:
            doc = make_valid_plan_doc([make_valid_stage_record(runtime_status=st)])
            validate_execution_plan(doc)

        # Database status
        for st in VALID_DATABASE_STATUS:
            doc = make_valid_plan_doc([make_valid_stage_record(database_status=st)])
            validate_execution_plan(doc)

        # Output status
        for st in VALID_OUTPUT_STATUS:
            doc = make_valid_plan_doc([make_valid_stage_record(output_status=st)])
            validate_execution_plan(doc)

        # Validation status
        for st in VALID_VALIDATION_STATUS:
            doc = make_valid_plan_doc([make_valid_stage_record(validation_status=st)])
            validate_execution_plan(doc)

    def test_missing_scalar_fields_serialize_to_null(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                checkpoint_reason=None,
                blocking_stage=None,
                reason_detail=None,
            )
        ])
        canonical = canonicalize_plan_data(doc)
        stg = canonical["plan"][0]
        self.assertIsNone(stg["checkpoint_reason"])
        self.assertIsNone(stg["blocking_stage"])
        self.assertIsNone(stg["reason_detail"])

        rendered = json.loads(render_json(doc))
        self.assertIsNone(rendered["plan"][0]["checkpoint_reason"])
        self.assertIsNone(rendered["plan"][0]["blocking_stage"])
        self.assertIsNone(rendered["plan"][0]["reason_detail"])

    def test_empty_lists_serialize_to_empty_list(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                invalidated_by=[],
                diagnostic_codes=[],
                dependencies=[],
                downstream=[],
            )
        ])
        canonical = canonicalize_plan_data(doc)
        stg = canonical["plan"][0]
        self.assertEqual(stg["invalidated_by"], [])
        self.assertEqual(stg["diagnostic_codes"], [])
        self.assertEqual(stg["dependencies"], [])
        self.assertEqual(stg["downstream"], [])

    def test_list_fields_ordered_and_duplicate_free(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                diagnostic_codes=["dup", "alpha", "dup", "beta"],
            )
        ])
        canonical = canonicalize_plan_data(doc)
        self.assertEqual(canonical["plan"][0]["diagnostic_codes"], ["dup", "alpha", "beta"])

    def test_cycle_detection_rejects_static_cycle(self):
        s1 = make_valid_stage_record("fastp", dependencies=["megahit"])
        s2 = make_valid_stage_record("megahit", dependencies=["fastp"])
        doc = make_valid_plan_doc([s1, s2])
        with self.assertRaises(ExecutionPlanInternalError):
            validate_execution_plan(doc)

    def test_missing_dependency_rejected(self):
        s1 = make_valid_stage_record("fastp", dependencies=["nonexistent_stage"])
        doc = make_valid_plan_doc([s1])
        with self.assertRaises(ExecutionPlanInternalError):
            validate_execution_plan(doc)

    def test_duplicated_stage_id_rejected(self):
        s1 = make_valid_stage_record("fastp")
        s2 = make_valid_stage_record("fastp")
        doc = make_valid_plan_doc([s1, s2])
        with self.assertRaises(ExecutionPlanInternalError):
            validate_execution_plan(doc)


class TestApplicabilityAndPrecedence(unittest.TestCase):
    """Test stage action determination, applicability, and precedence rules."""

    def test_disabled_stage_resolves_to_disabled(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "fastp",
                configured_enabled=False,
                planned_action="disabled",
                reason_code="stage-disabled",
                required=False,
            )
        ])
        validate_execution_plan(doc)
        canonical = canonicalize_plan_data(doc)
        self.assertEqual(canonical["plan"][0]["planned_action"], "disabled")
        self.assertFalse(canonical["plan"][0]["required"])

    def test_not_applicable_stage_resolves_to_not_applicable(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "reassembly",
                configured_enabled=True,
                applicability="not-applicable",
                planned_action="not-applicable",
                reason_code="condition-not-met",
                required=False,
            )
        ])
        validate_execution_plan(doc)
        canonical = canonicalize_plan_data(doc)
        self.assertEqual(canonical["plan"][0]["planned_action"], "not-applicable")
        self.assertFalse(canonical["plan"][0]["required"])

    def test_precedence_disabled_and_not_applicable_precede_resource_inspection(self):
        # Even if database is missing and runtime is unavailable, a disabled or not-applicable
        # stage stays disabled or not-applicable, NOT blocked.
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "gtdbtk",
                configured_enabled=False,
                applicability="applicable",
                database_status="not-configured",
                runtime_status="unavailable",
                planned_action="disabled",
                reason_code="stage-disabled",
                required=False,
            ),
            make_valid_stage_record(
                "reassembly",
                configured_enabled=True,
                applicability="not-applicable",
                database_status="not-configured",
                runtime_status="unavailable",
                planned_action="not-applicable",
                reason_code="condition-not-met",
                required=False,
            ),
        ])
        validate_execution_plan(doc)
        self.assertTrue(is_plan_executable(doc["plan"]))

    def test_indeterminate_applicability_blocks_execution(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "binning",
                applicability="indeterminate",
                planned_action="blocked",
                reason_code="indeterminate-applicability",
                required=True,
            )
        ])
        validate_execution_plan(doc)
        self.assertFalse(is_plan_executable(doc["plan"]))


class TestCheckpointsAndValidation(unittest.TestCase):
    """Test checkpoint handling, validation states, and resource blockages."""

    def test_valid_checkpoint_resolves_to_reuse(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "fastp",
                checkpoint_status="valid",
                output_status="available",
                validation_status="valid",
                planned_action="reuse",
                reason_code="checkpoint-valid",
            )
        ])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["reuse"], 1)
        self.assertEqual(summary["run"], 0)

    def test_missing_checkpoint_resolves_to_run(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "fastp",
                checkpoint_status="missing",
                planned_action="run",
                reason_code="checkpoint-missing",
            )
        ])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["run"], 1)

    def test_missing_config_manifest_with_present_done_resolves_to_run(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "fastp",
                checkpoint_status="legacy",
                checkpoint_reason="manifest-missing",
                planned_action="run",
                reason_code="checkpoint-legacy",
            )
        ])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["run"], 1)

    def test_corrupted_checkpoint_resolves_to_run(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "fastp",
                checkpoint_status="invalid",
                checkpoint_reason="manifest-corrupted",
                planned_action="run",
                reason_code="checkpoint-invalid",
            )
        ])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["run"], 1)

    def test_validation_failure_on_reusable_checkpoint_resolves_to_run(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "fastp",
                checkpoint_status="valid",
                validation_status="invalid",
                planned_action="run",
                reason_code="validation-failed",
            )
        ])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["run"], 1)

    def test_unusable_required_runtime_blocks_stage(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "megahit",
                runtime_status="unavailable",
                planned_action="blocked",
                reason_code="runtime-unavailable",
                required=True,
            )
        ])
        validate_execution_plan(doc)
        self.assertFalse(is_plan_executable(doc["plan"]))

    def test_missing_required_database_blocks_stage(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "kraken2",
                database_status="not-configured",
                planned_action="blocked",
                reason_code="database-not-configured",
                required=True,
            )
        ])
        validate_execution_plan(doc)
        self.assertFalse(is_plan_executable(doc["plan"]))

    def test_missing_input_file_blocks_stage(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "fastp",
                input_status="unavailable",
                planned_action="blocked",
                reason_code="input-missing",
                required=True,
            )
        ])
        validate_execution_plan(doc)
        self.assertFalse(is_plan_executable(doc["plan"]))


class TestForceAndRestart(unittest.TestCase):
    """Test requested_action preservation and cascade rules for force and restart."""

    def test_requested_action_preserved(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record("fastp", requested_action="force", planned_action="force"),
            make_valid_stage_record("megahit", requested_action="restart", planned_action="restart"),
            make_valid_stage_record("mapping", requested_action="normal", planned_action="run"),
        ])
        validate_execution_plan(doc)
        self.assertEqual(doc["plan"][0]["requested_action"], "force")
        self.assertEqual(doc["plan"][1]["requested_action"], "restart")
        self.assertEqual(doc["plan"][2]["requested_action"], "normal")

    def test_force_stage_target_and_downstream(self):
        s_target = make_valid_stage_record(
            "fastp",
            requested_action="force",
            planned_action="force",
            reason_code="forced-by-user",
            downstream=["megahit"],
        )
        s_downstream = make_valid_stage_record(
            "megahit",
            requested_action="normal",
            planned_action="run",
            reason_code="upstream-forced",
            dependencies=["fastp"],
            diagnostic_codes=["upstream-forced:fastp"],
            invalidated_by=["fastp"],
        )
        doc = make_valid_plan_doc([s_target, s_downstream])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["force"], 1)
        self.assertEqual(summary["run"], 1)

    def test_restart_from_boundary_and_downstream(self):
        s_target = make_valid_stage_record(
            "megahit",
            requested_action="restart",
            planned_action="restart",
            reason_code="restart-boundary",
            downstream=["mapping"],
        )
        s_downstream = make_valid_stage_record(
            "mapping",
            requested_action="normal",
            planned_action="run",
            reason_code="downstream-of-restart",
            dependencies=["megahit"],
            diagnostic_codes=["downstream-of-restart:megahit"],
            invalidated_by=["megahit"],
        )
        doc = make_valid_plan_doc([s_target, s_downstream])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["restart"], 1)
        self.assertEqual(summary["run"], 1)

    def test_force_or_restart_on_blocked_stage_remains_blocked(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record(
                "kraken2",
                requested_action="force",
                database_status="not-configured",
                planned_action="blocked",
                reason_code="database-not-configured",
                diagnostic_codes=["forced-by-user"],
            ),
            make_valid_stage_record(
                "checkm2",
                requested_action="restart",
                database_status="not-configured",
                planned_action="blocked",
                reason_code="database-not-configured",
                diagnostic_codes=["restart-target"],
            ),
        ])
        validate_execution_plan(doc)
        for stg in doc["plan"]:
            self.assertEqual(stg["planned_action"], "blocked")


class TestGraphAndInvalidation(unittest.TestCase):
    """Test branch isolation and cascade invalidation."""

    def test_kraken2_invalidation_does_not_invalidate_mag_branch(self):
        # Kraken2 downstream is empty or only taxonomy-summary; it must not invalidate megahit or binning
        s_kraken = make_valid_stage_record("kraken2", planned_action="force", requested_action="force", downstream=[])
        s_megahit = make_valid_stage_record("megahit", checkpoint_status="valid", planned_action="reuse")
        s_binning = make_valid_stage_record("binning", checkpoint_status="valid", planned_action="reuse", dependencies=["megahit"])
        doc = make_valid_plan_doc([s_kraken, s_megahit, s_binning])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["reuse"], 2)
        self.assertEqual(summary["force"], 1)

    def test_checkm2_invalidation_invalidates_filtering_and_downstream_mag_steps(self):
        s_checkm2 = make_valid_stage_record(
            "checkm2", planned_action="force", requested_action="force", downstream=["mag_filtering"]
        )
        s_filter = make_valid_stage_record(
            "mag_filtering",
            planned_action="run",
            reason_code="upstream-forced",
            dependencies=["checkm2"],
            downstream=["gtdbtk"],
            invalidated_by=["checkm2"],
        )
        s_gtdbtk = make_valid_stage_record(
            "gtdbtk",
            planned_action="run",
            reason_code="upstream-forced",
            dependencies=["mag_filtering"],
            invalidated_by=["checkm2"],
        )
        doc = make_valid_plan_doc([s_checkm2, s_filter, s_gtdbtk])
        validate_execution_plan(doc)
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(summary["force"], 1)
        self.assertEqual(summary["run"], 2)


class TestExecutabilityAndSummary(unittest.TestCase):
    """Test plan-level executability and summary counts."""

    def test_executable_true_when_all_required_stages_executable(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record("fastp", planned_action="run", required=True),
            make_valid_stage_record("megahit", planned_action="reuse", required=True),
            make_valid_stage_record("reassembly", planned_action="disabled", required=False),
        ])
        self.assertTrue(is_plan_executable(doc["plan"]))

    def test_executable_false_when_required_stage_blocked(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record("fastp", planned_action="run", required=True),
            make_valid_stage_record("megahit", planned_action="blocked", required=True),
        ])
        self.assertFalse(is_plan_executable(doc["plan"]))

    def test_non_required_blocked_stage_does_not_make_plan_non_executable(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record("fastp", planned_action="run", required=True),
            make_valid_stage_record("kraken2", planned_action="blocked", required=False),
        ])
        self.assertTrue(is_plan_executable(doc["plan"]))

    def test_summary_counts_match_individual_stage_actions(self):
        doc = make_valid_plan_doc([
            make_valid_stage_record("fastp", planned_action="run"),
            make_valid_stage_record("megahit", planned_action="reuse", checkpoint_status="valid"),
            make_valid_stage_record("kraken2", planned_action="force", requested_action="force"),
            make_valid_stage_record("mapping", planned_action="restart", requested_action="restart"),
            make_valid_stage_record("binning", planned_action="disabled", configured_enabled=False, required=False),
            make_valid_stage_record("reassembly", planned_action="not-applicable", applicability="not-applicable", required=False),
            make_valid_stage_record("checkm2", planned_action="blocked", required=True),
        ])
        summary = compute_plan_summary(doc["plan"])
        self.assertEqual(sum(summary.values()), 7)
        self.assertEqual(summary["run"], 1)
        self.assertEqual(summary["reuse"], 1)
        self.assertEqual(summary["force"], 1)
        self.assertEqual(summary["restart"], 1)
        self.assertEqual(summary["disabled"], 1)
        self.assertEqual(summary["not_applicable"], 1)
        self.assertEqual(summary["blocked"], 1)


class TestCLISideEffectsAndInvocations(unittest.TestCase):
    """Integration tests running run_metagenomics_pipeline.sh CLI."""

    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="metanexus_test_")
        self.outdir = os.path.join(self.test_dir, "test_outdir")

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run_pipeline(self, args, check=False):
        cmd = [str(SCRIPT_PATH)] + args
        env = dict(os.environ)
        # Avoid SIF container resolution issues in testing environment
        res = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            env=env,
        )
        if check and res.returncode != 0:
            raise subprocess.CalledProcessError(res.returncode, cmd, res.stdout, res.stderr)
        return res

    def test_plan_does_not_create_outdir(self):
        self.assertFalse(os.path.exists(self.outdir))
        res = self.run_pipeline([
            "--plan",
            "-o", self.outdir,
            "-i", str(REPO_ROOT / "config"),  # dummy input dir
        ])
        self.assertFalse(
            os.path.exists(self.outdir),
            f"OUTDIR {self.outdir} was created during --plan! Contract violation."
        )

    def test_plan_does_not_write_checkpoints_or_manifests(self):
        res = self.run_pipeline([
            "--plan",
            "-o", self.outdir,
        ])
        # OUTDIR must not even exist
        self.assertFalse(os.path.exists(self.outdir))
        # Check no persistent files leaked into repo root
        for leaked in ["resume_plan.tsv", "stage_status.tsv", "run_overview.tsv", "errors.tsv", "last_error.json"]:
            self.assertFalse(
                (REPO_ROOT / leaked).exists(),
                f"Leaked persistent file {leaked} in repo root during --plan"
            )

    def test_plan_options_require_plan(self):
        # --plan-format requires --plan
        res = self.run_pipeline(["--plan-format", "json"])
        self.assertEqual(res.returncode, 2)
        self.assertIn("--plan", res.stderr)

        # --plan-details requires --plan
        res = self.run_pipeline(["--plan-details", "summary"])
        self.assertEqual(res.returncode, 2)
        self.assertIn("--plan", res.stderr)

        # --strict-plan requires --plan
        res = self.run_pipeline(["--strict-plan"])
        self.assertEqual(res.returncode, 2)
        self.assertIn("--plan", res.stderr)

    def test_plan_help_exits_zero(self):
        res = self.run_pipeline(["--plan", "--help"])
        self.assertEqual(res.returncode, 0)
        self.assertIn("Execution Plan Inspection", res.stdout)
        self.assertIn("--plan-format", res.stdout)
        self.assertIn("--plan-details", res.stdout)
        self.assertIn("--strict-plan", res.stdout)

    def test_plan_incompatible_modes(self):
        # --plan + --dry-run
        res = self.run_pipeline(["--plan", "--dry-run"])
        self.assertEqual(res.returncode, 2)
        self.assertIn("Cannot combine --plan with", res.stderr)

        # --plan + --version
        res = self.run_pipeline(["--plan", "--version"])
        self.assertEqual(res.returncode, 2)
        self.assertIn("Cannot combine --plan with", res.stderr)

        # --plan + --versions
        res = self.run_pipeline(["--plan", "--versions"])
        self.assertEqual(res.returncode, 2)
        self.assertIn("Cannot combine --plan with", res.stderr)

    def test_plan_json_output(self):
        res = self.run_pipeline([
            "--plan",
            "--plan-format", "json",
            "-o", self.outdir,
        ])
        self.assertIn(res.returncode, [0, 1])
        try:
            data = json.loads(res.stdout)
        except Exception as e:
            self.fail(f"Output was not valid JSON: {e}\nOutput was:\n{res.stdout}\nStderr:\n{res.stderr}")
        self.assertEqual(data.get("schema_version"), "1")
        self.assertEqual(data.get("pipeline", {}).get("name"), PIPELINE_NAME)
        self.assertEqual(data.get("pipeline", {}).get("version"), PIPELINE_VERSION)
        self.assertEqual(data.get("pipeline", {}).get("release_status"), PIPELINE_RELEASE_STATUS)
        self.assertIsInstance(data.get("plan"), list)

    def test_plan_exit_code_zero_when_blocked(self):
        # Default run without inputs has blocked stages, but --plan without --strict-plan exits 0
        res = self.run_pipeline([
            "--plan",
            "-o", self.outdir,
        ])
        self.assertEqual(res.returncode, 0, f"Expected 0, got {res.returncode}. Stderr: {res.stderr}")

    def test_strict_plan_exit_code_one_when_not_executable(self):
        # Default run without inputs has blocked required stages, so executable=false -> exit 1
        res = self.run_pipeline([
            "--plan",
            "--strict-plan",
            "-o", self.outdir,
        ])
        self.assertEqual(
            res.returncode, 1,
            f"Expected exit code 1 for blocked plan under --strict-plan, got {res.returncode}. Stderr: {res.stderr}"
        )


if __name__ == "__main__":
    unittest.main()
