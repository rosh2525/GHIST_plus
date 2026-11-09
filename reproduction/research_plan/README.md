# GHIST+ research implementation plan

Current revision: 31 October 2026. This folder contains planning documents, not an implemented or trained new model.

- [Complete implementation plan](GHIST_PLUS_IMPLEMENTATION_PLAN.md): audited baseline, data, selective interval-supervision candidate, correctness checks, evaluation, gates, resources and schedule.
- [Dataset registry](dataset_sources.json): source accessions and provisional donor roles; not a download or image-QC receipt.
- [Experiment matrix](EXPERIMENT_MATRIX.csv): baseline, target-cleaning, masking and proposed-method comparisons.
- [Execution specification](execution_plan.json): dependencies, gates and budget assumptions; not a runnable training config.
- [Previous correspondence-first plan](archive/GHIST_PLUS_IMPLEMENTATION_PLAN_2026-10-30.md): preserved for provenance.

The main candidate is conditional on the two-week audit. Corrected RNA and nuclear-compartment counts are not treated as clean whole-cell truth. The new objective must beat ordinary target cleaning and masking on independent evidence. Exact manuscript reproduction remains separate in [the existing report](../REPORT.md).

- [Gemini review decisions](REVIEW_GEMINI_FLASH_2026-10-31.md): accepted safeguards, corrected mathematical claims, and analytic checks.
