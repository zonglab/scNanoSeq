# Replace an existing local clone with this refreshed working tree

The distributed ZIP intentionally excludes `.git`. This lets you keep the Git history/remotes from your existing clone while replacing only the working-tree files.

## Safest method: preserve `.git`, replace everything else

Assume:

```text
/path/to/scNanoSeq/            # your existing git clone
/path/to/scNanoSeq_refreshed/  # the extracted refreshed folder
```

First inspect your existing clone:

```bash
cd /path/to/scNanoSeq
git status
```

If you have local work you care about, commit it, stash it, or make a backup before continuing.

Then replace the working tree while preserving `.git`:

```bash
rsync -av --delete --exclude='.git/' /path/to/scNanoSeq_refreshed/ /path/to/scNanoSeq/
```

Now inspect exactly what changed:

```bash
cd /path/to/scNanoSeq
git status
git diff -- README.md
git diff -- 1_Bulk_WGS 2_Somatic_mutation_calling 2_python_script 3_Phylogeny 4_CNV
```

Because the uploaded snapshot had already reorganized the historical directory names into numbered directories, `git status` will show many old paths as deleted and numbered paths as new until you commit. Git may later recognize some of these as renames automatically.

Run the repository checks:

```bash
bash scripts/repo_check.sh
```

If `samtools` and `pysam` are installed, run the core synthetic demo:

```bash
bash demo/run_core_demo.sh
```

After reviewing the changes:

```bash
git add -A
git status
git diff --cached
```

Commit only when the staged diff matches what you intend to publish:

```bash
git commit -m "Document scNanoSeq pipeline and add minimal synthetic demo"
git push origin main
```

## If you want to keep the old directory layout instead

Do not use `rsync --delete`. Copy only the new documentation/demo files and manually port the two small script fixes. The refreshed package follows the numbered layout present in the uploaded ZIP, not the historical layout stored at Git commit `2f1c8cb`.

## Files that should not be overwritten blindly

Before publishing, pay special attention to:

- any locally recovered `gtbulksearch_py3.py`;
- any real project-specific path configuration;
- an existing top-level `LICENSE` or `CITATION.cff` that was not present in the uploaded ZIP; and
- any unpublished analysis scripts or data that exist only in your local clone.
