# Safely applying this review-ready folder to an existing Git clone

The review-ready bundle intentionally does **not** contain `.git`. Keep the `.git` directory from your existing clone so your history/remotes remain intact.

## Recommended Linux/macOS/HPC procedure

Assume:

- your existing clone is `/path/to/scNanoSeq`;
- this bundle has been unzipped as `/path/to/scNanoSeq_review_ready`.

First protect your current work:

```bash
cd /path/to/scNanoSeq
git status
# If you have uncommitted work, commit it or stash it before continuing.
git switch -c nature-software-checklist
cd ..
cp -a scNanoSeq scNanoSeq.backup-before-checklist
```

Then mirror the review-ready contents **while preserving `.git`**:

```bash
rsync -av --delete --exclude='.git' \
  /path/to/scNanoSeq_review_ready/ \
  /path/to/scNanoSeq/
```

`--delete` makes the working tree match the supplied bundle exactly. That is why the backup step exists. After copying:

```bash
cd /path/to/scNanoSeq
git status
git add -A
git diff --cached --stat
git diff --cached
```

Review the diff before committing. In particular, confirm whether removal of the historical `Simulation/` and `VAF_fingerprint/` directories is intentional. Those directories exist on the current public `main` branch but were absent from the uploaded working tree used to prepare this bundle.

When satisfied:

```bash
git commit -m "Prepare scNanoSeq for software review"
git push -u origin nature-software-checklist
```

Merge to `main` only after you have supplied the remaining AUTHOR ACTION REQUIRED items.

## Safer non-destructive copy

If you do not want any existing file deleted automatically, omit `--delete`:

```bash
rsync -av --exclude='.git' /path/to/scNanoSeq_review_ready/ /path/to/scNanoSeq/
```

This will leave obsolete directories/files in place, so you must remove or reconcile them manually before committing.

## Windows PowerShell / Robocopy

From a PowerShell prompt, after creating a backup and a new Git branch:

```powershell
robocopy C:\path\to\scNanoSeq_review_ready C:\path\to\scNanoSeq /MIR /XD .git
```

`/MIR` deletes destination files not present in the source bundle, so use it only after making the backup and checking that you want an exact mirror.

## Validation after replacement

Run:

```bash
python3 demo/run_demo.py
python3 scripts/reviewer_preflight.py
```

The demo should pass. The preflight is expected to remain non-zero until the project license and missing `gtbulksearch_py3.py` are supplied; those failures are deliberate rather than hidden.
