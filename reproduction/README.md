# GHIST+ manuscript-only reproduction

Project scope: **GHIST+**, the latest direct successor to GHIST verified on **22 October 2026**. Reproduce the public manuscript before attempting an improvement.

- [Current report](REPORT.md): paper-reported values, measured release results, and unresolved differences.
- [Step-by-step workflow](WORKFLOW.md): datasets, preprocessing, architecture, experiments, settings and commands.
- [Paper protocol](paper_protocol/protocol.json): metric whitelist, exact PDF version/hash and experiment definitions.
- [Paper metric values](paper_protocol/paper_values.csv): numbers transcribed from the manuscript with provenance.

**The released notebooks ran successfully, but exact numerical reproduction of the paper has not been achieved.** The public release and manuscript differ in values, normalization, folds and figure contents. The current report preserves those distinctions.

Run from the `GHIST_plus` repository root using the existing isolated environment:

```sh
../work/ghist_plus_reproduction/.venv/bin/python reproduction/fetch_paper.py
../work/ghist_plus_reproduction/.venv/bin/python reproduction/prepare_paper_templates.py
../work/ghist_plus_reproduction/.venv/bin/python reproduction/paper_only.py
```

The earlier broad notebook review is [archived](archive/2026-10-15_released_notebook_review.md). Its Phoenix, ATERA/Breast5K and additional cohort comparisons are outside the current manuscript-only scope. Existing downloads and original outputs were preserved.

The official code and notebook cells are unchanged. The new templates apply settings explicitly stated in the manuscript and flag missing original fold assignments. They are not represented as archived paper configs. Fresh training requires prepared real images/RNA, approved UNI2-H access and a suitable model environment.

Source code remains pinned at `c5c895174c19903cc27d2661934b9b726bdb9306`; data revision `37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de`. The local branch is `reproduction/ghist-plus-local`; no changes are pushed upstream. See the original README and third-party notices for applicable code/data/encoder terms.
