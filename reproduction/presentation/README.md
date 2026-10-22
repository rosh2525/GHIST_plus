# GHIST+ workflow presentation

Six editable slides, prepared 22 October 2026. `GHIST_plus_workflow.pptx` is the reviewed deck. `presenter_notes.md` contains the detailed slide explanations, source links and conditional training commands; the same notes are embedded in the PPTX.

## Scope

The completed work is evaluation of released predictions against real measured RNA, plus verification of recorded notebook execution. No fresh model training or image inference is claimed. Exact manuscript values and the complete original five-fold, pancreas and GTEx experiments remain incomplete.

Slides cover status, the actual dataset cohorts, the model workflow, local execution, paper versus released-code PCC, and the remaining reproduction workflow. The comparison keeps the differing values visible. Do not call PCC percent accuracy.

## Screenshots and data provenance

`evidence/*.jpg` are real browser captures of local pages displaying actual command stdout. They are not native Terminal screenshots. `transcripts.json` records each command, stdout/stderr, return code, capture time and the evidence-script SHA256. The 15 October notebook results are recorded execution results; the provenance/cohort/PCC check was rerun on 22 October. No shell history, credentials or private terminal sessions are included.

From the repository root, with the existing analysis environment:

```sh
../work/ghist_plus_reproduction/.venv/bin/python reproduction/presentation_evidence.py data
../work/ghist_plus_reproduction/.venv/bin/python reproduction/presentation_evidence.py checks
../work/ghist_plus_reproduction/.venv/bin/python reproduction/presentation_evidence.py compare
../work/ghist_plus_reproduction/.venv/bin/python reproduction/presentation_evidence.py status
```

`source_tables/` freezes the counts and metric comparison used in this deck. The input RNA and released prediction arrays remain in the ignored data folder. No synthetic data is included.

`build.mjs` is the Artifact Tool source used for this local deck. It requires the Codex bundled Node runtime, `@oai/artifact-tool`, the presentation skill finalizer, and the local workspace paths recorded at its top. It reads `work/ghist_workflow_deck/evidence` and the paper-only result tables. Choose new candidate/final/receipt filenames before rebuilding: the finalizer deliberately does not overwrite a previous final artifact. `capture_evidence.py` records real stdout into HTML/transcripts; the JPEGs were then captured with the browser tool.

## References

- Paper v1: https://doi.org/10.64898/2026.06.09.731252
- Official code: https://github.com/SydneyBioX/GHIST_plus
- Processed real data: https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle
- Training workflow and caveats: ../WORKFLOW.md

Code pin: c5c895174c19903cc27d2661934b9b726bdb9306. Dataset pin: 37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de.
