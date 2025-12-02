# Validation Engine Roadmap

## Overview

This document outlines the plan to integrate and complete the module artifact validation feature for omnibenchmark. The validation engine validates outputs from benchmark modules to ensure they meet quality and format requirements.

## Current State

### What Exists

1. **New Validation Engine** (`omnibenchmark/validators/`)
   - Located in: `omnibenchmark/validators/engine.py`
   - Uses `returns` library for robust error handling
   - Implements validation rules: `is_file`, `is_dir`, `matches`, `checksum`, `has_column`, `not_na`, `has_key`
   - Returns `ValidationResult` and `ModuleValidationResult` dataclasses
   - Tests in `tests/validators/test_engine.py`

2. **Old Artifact Validation** (`omnibenchmark/artifact_validation/`)
   - More comprehensive with extensible validator/loader architecture
   - Uses Pydantic models for validation results
   - Has `ValidationStatus` enum with severity levels: PASSED, FAILED, SKIPPED, ERROR
   - Includes data loaders for CSV, JSON, etc.
   - Has `ValidationReport` for benchmark-wide results
   - **Not committed to the repository**

3. **Design Documentation**
   - `docs/design/002-module-artifact-validation.md`
   - Outlines the validation.yaml specification
   - Describes stage-based validation approach

## Goals

### Phase 1: Severity Levels and Warnings (Priority: HIGH)

**Goal**: Add `--warn` flag to convert strict validation failures into warnings

**Tasks**:

1. **Add Severity Levels to Validation Rules**
   - [ ] Define severity enum: `ERROR`, `WARN`, `INFO`
   - [ ] Extend `validation.yaml` spec to support severity per rule:
     ```yaml
     stages:
       dataset:
         rules:
           - file_pattern: "*.csv"
             validations:
               - type: not_empty
                 severity: error  # default
               - type: has_columns
                 columns: ["id", "value"]
                 severity: warn
     ```
   - [ ] Update validation rule types to include optional `severity` field

2. **Update ValidationResult to Track Severity**
   - [ ] Add `severity` field to `ValidationResult` dataclass
   - [ ] Distinguish between errors (block execution) and warnings (log only)

3. **Add CLI Flag for Strict Mode**
   - [ ] Add `--warn` / `--strict` flag to relevant CLI commands
   - [ ] Default: strict mode (all failures are errors)
   - [ ] With `--warn`: convert ERROR → WARN (except critical validations)
   - [ ] Implement in workflow execution to respect this flag

4. **Reporting**
   - [ ] Update result formatting to show severity levels
   - [ ] Color-code output: RED for errors, YELLOW for warnings, GREEN for passed
   - [ ] Summary should count errors vs warnings separately

**Files to Modify**:
- `omnibenchmark/validators/engine.py`
- `omnibenchmark/cli/run.py` (or relevant CLI module)
- `omnibenchmark/workflow/snakemake/snakemake.py`

---

### Phase 2: Per-Module validation.yaml Support (Priority: HIGH)

**Goal**: Enable modules to define their own validation rules

**Tasks**:

1. **Module validation.yaml Discovery**
   - [ ] Define standard location: `<module_repo>/validation.yaml`
   - [ ] Add logic to clone/fetch module repo and look for validation.yaml
   - [ ] If not found, use benchmark-level validation rules (if any)

2. **Validation Rule Resolution**
   - [ ] Implement rule merging strategy:
     - Module-level rules override benchmark-level rules
     - Or: Both apply (benchmark + module)
   - [ ] Document the precedence/merging behavior

3. **Integration with Module Execution**
   - [ ] After module execution, run validation before marking output as complete
   - [ ] Hook into `run_module.py` or workflow execution
   - [ ] If validation fails in strict mode, mark module execution as failed

4. **Validation Rule Templates**
   - [ ] Create example validation.yaml for common module types:
     - Dataset modules
     - Method modules  
     - Metric modules
   - [ ] Add to omnibenchmark documentation

**Files to Modify**:
- `omnibenchmark/workflow/snakemake/scripts/run_module.py`
- `omnibenchmark/io/code.py` (module cloning logic)
- `omnibenchmark/validators/engine.py`

---

### Phase 3: Consolidate Old and New Validation Code (Priority: MEDIUM)

**Goal**: Merge useful patterns from old artifact_validation into new validators

**Tasks**:

1. **Review Old vs New Architecture**
   - [ ] Compare `omnibenchmark/artifact_validation/` with `omnibenchmark/validators/`
   - [ ] Identify valuable patterns from old code:
     - Extensible validator registry
     - Data loader abstraction
     - Comprehensive status enum
     - Pydantic models vs dataclasses

2. **Migrate Valuable Components**
   - [ ] Consider adding data loader abstraction to new engine
   - [ ] Evaluate if Pydantic models offer benefits over dataclasses
   - [ ] Port any useful validators from old code
   - [ ] Keep or adapt `ValidationStatus` enum (PASSED/FAILED/SKIPPED/ERROR)

3. **Deprecate or Remove Old Code**
   - [ ] Once migration complete, remove `omnibenchmark/artifact_validation/`
   - [ ] Update any references to old validation system
   - [ ] Add deprecation warnings if needed for gradual migration

**Files to Review**:
- `omnibenchmark/artifact_validation/base.py`
- `omnibenchmark/artifact_validation/engine.py`
- `omnibenchmark/artifact_validation/loaders/`
- `omnibenchmark/artifact_validation/validators/`

---

### Phase 4: Workflow Integration (Priority: DEFERRED)

**Goal**: Automatically validate module outputs during workflow execution

> **Note**: This phase is DEFERRED for now. Focus is on standalone validation tools first. However, we keep these notes for future consideration.

**Future Considerations**:

1. **Post-Execution Validation Hook**
   - [ ] Add validation step after module execution in `run_module.py`
   - [ ] Load validation rules from module or benchmark
   - [ ] Run validation on module outputs
   - [ ] Log results to validation report
   - **Feasibility**: ✅ Straightforward to add after module execution completes
   - **Challenge**: Need to handle validation failures without breaking workflow

2. **Continue-on-Error Interaction**
   - [ ] When `--continue-on-error` is set:
     - Still run validation
     - Log failures but don't stop workflow
   - [ ] When strict mode (default):
     - Validation failures cause module to fail
     - Don't touch empty files (current bug)
   - **Feasibility**: ✅ Can integrate with existing error handling
   - **Challenge**: Distinguish validation failures from execution failures

3. **Validation Report Collection**
   - [ ] Collect validation results from all modules
   - [ ] Generate benchmark-wide validation report
   - [ ] Store report in output directory (e.g., `validation_report.json`)
   - **Feasibility**: ✅ Can aggregate results during workflow execution
   - **Challenge**: Need to persist state across Snakemake rules

4. **Snakemake Integration Options**
   - Option A: Integrate directly into module execution rules
     - **Pros**: Single step, simpler workflow
     - **Cons**: Validation errors harder to distinguish from execution errors
   - Option B: Add validation as separate Snakemake rules
     - **Pros**: Clear separation, can skip validation if needed
     - **Cons**: More complex DAG, harder to make validation optional
   - Option C: Post-workflow validation pass
     - **Pros**: Non-invasive, works with existing workflows
     - **Cons**: Doesn't prevent metric collectors from using invalid data
   - **Recommendation**: Start with Option C (post-workflow), consider Option A later

5. **Integration Points**
   - `run_module.py`: After `execution()` call, before `sys.exit(0)`
   - `rule_node.smk`: Could add `validation.yaml` as input dependency
   - Snakemake config: Add `validation_mode` option (strict/warn/skip)
   - **Feasibility**: ✅ All integration points are accessible

**Why Deferred**:
- Standalone validation tools are more immediately useful
- Allows module authors to test validation rules independently
- Workflow integration needs more design discussion
- Risk of breaking existing benchmarks during development

**Files That Would Be Modified** (for future reference):
- `omnibenchmark/workflow/snakemake/scripts/run_module.py`
- `omnibenchmark/workflow/snakemake/rules/rule_node.smk`
- `omnibenchmark/workflow/snakemake/snakemake.py`

---

### Phase 5: Additional Validation Rules (Priority: LOW)

**Goal**: Expand validation capabilities

**Potential Validators**:

1. **Statistical Validators**
   - [ ] `value_range`: Check numeric values are within range
   - [ ] `no_duplicates`: Ensure unique values in column
   - [ ] `distribution_check`: Basic statistical checks

2. **Format Validators**
   - [ ] `valid_json_schema`: Validate against JSON schema
   - [ ] `valid_format`: Check file format (e.g., gzipped correctly)
   - [ ] `encoding_check`: Verify file encoding

3. **Relational Validators**
   - [ ] `references_exist`: Check foreign key relationships
   - [ ] `consistent_dimensions`: Verify matrix/array dimensions match

**Files to Modify**:
- `omnibenchmark/validators/engine.py`
- Add new handler functions for each validator

---

### Phase 6: CLI and User Interface (Priority: MEDIUM)

**Goal**: Provide user-friendly validation tools

**Tasks**:

1. **Standalone Validation Command**
   - [ ] Add `ob validate` command
   - [ ] Usage: `ob validate <output_dir> --config validation.yaml`
   - [ ] Can be run independently of workflow execution

2. **Validation Status Command**
   - [ ] Add to `ob status` command to show validation status
   - [ ] Display which modules passed/failed validation

3. **Validation Report Command**
   - [ ] Add `ob report validation` subcommand
   - [ ] Display validation report with formatting
   - [ ] Support multiple output formats (text, JSON, HTML)

**Files to Modify**:
- `omnibenchmark/cli/validate.py` (new file)
- `omnibenchmark/cli/status.py`
- `omnibenchmark/cli/report.py`

---

## Implementation Order

### Sprint 1: Foundation (Weeks 1-2)
1. Phase 1: Severity levels and --warn flag
2. Phase 6 (partial): Basic `ob validate` CLI command

### Sprint 2: Module Support (Weeks 3-4)
1. Phase 2: Per-module validation.yaml support
2. Phase 3: Consolidate old and new code

### Sprint 3: Polish (Weeks 5-6)
1. Phase 6 (complete): Full CLI and user interface
2. Documentation and examples
3. Integration tests for validation

### Future Work
- Phase 4: Workflow integration (deferred, see notes in Phase 4)
- Phase 5: Additional validation rules (as needed)

---

## Open Questions

1. **Severity Behavior**:
   - Should `--warn` convert ALL errors to warnings, or only specific rules?
   - Should some rules always be errors (e.g., file not found)?

2. **Module vs Benchmark Rules**:
   - If both exist, which takes precedence?
   - Should we merge them or let module rules completely override?

3. **Validation Timing**:
   - Validate immediately after module execution?
   - Or as a separate Snakemake step?
   - What about metric collectors - do they need validation too?

4. **Empty File Handling**:
   - Current bug: touched files from failed modules are empty
   - Should validation detect and mark these as invalid?
   - Or should we fix the touching behavior instead?

5. **Backward Compatibility**:
   - Do we need to support old benchmarks without validation?
   - Should validation be opt-in or opt-out?

---

## Success Criteria

✅ Phase 1 Complete When:
- `--warn` flag works in CLI
- Validation results show severity levels
- Errors vs warnings are clearly distinguished

✅ Phase 2 Complete When:
- Modules can define validation.yaml
- Rules are discovered and applied automatically
- Validation runs during workflow execution

✅ Phase 4 Complete When:
- All module outputs are validated automatically
- Validation reports are generated per benchmark run
- Metric collectors don't process invalid outputs

✅ Overall Feature Complete When:
- Integration tests pass with validation enabled
- Documentation is complete with examples
- e2e tests include validation checks
- Old artifact_validation code is deprecated/removed

---

## Notes

- The new validation engine (`omnibenchmark/validators/`) is simpler and uses functional error handling
- The old validation engine (`omnibenchmark/artifact_validation/`) is more extensible but not committed
- Consider best of both approaches when consolidating
- Severity levels are critical for developer experience (don't block on minor issues)
- Per-module validation.yaml is essential for module authors to define their contracts
