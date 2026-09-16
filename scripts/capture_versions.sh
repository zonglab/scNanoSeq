#!/usr/bin/env bash
set -u

printf 'scNanoSeq software-version capture\n'
printf 'date\t%s\n' "$(date -Iseconds 2>/dev/null || date)"
printf 'git_commit\t%s\n' "$(git rev-parse HEAD 2>/dev/null || echo NA)"
printf 'python3\t%s\n' "$(python3 --version 2>&1 || echo NA)"

for tool in cutadapt bwa samtools bcftools bedtools seqtk picard java Rscript cellphy.sh; do
    if command -v "$tool" >/dev/null 2>&1; then
        case "$tool" in
            bwa) out=$(bwa 2>&1 | head -n 3 | tr '\n' ' ') ;;
            bedtools) out=$(bedtools --version 2>&1 | head -n 1) ;;
            seqtk) out=$(seqtk 2>&1 | head -n 1) ;;
            picard) out=$(picard -h 2>&1 | head -n 2 | tr '\n' ' ') ;;
            java) out=$(java -version 2>&1 | head -n 1) ;;
            Rscript) out=$(Rscript --version 2>&1 | head -n 1) ;;
            cellphy.sh) out=$(cellphy.sh 2>&1 | head -n 2 | tr '\n' ' ') ;;
            *) out=$($tool --version 2>&1 | head -n 1) ;;
        esac
        printf '%s\t%s\n' "$tool" "$out"
    else
        printf '%s\tNOT_FOUND\n' "$tool"
    fi
done

python3 - <<'PY'
for package in ("numpy", "pandas", "pysam"):
    try:
        mod = __import__(package)
        print(f"python:{package}\t{getattr(mod, '__version__', 'unknown')}")
    except Exception:
        print(f"python:{package}\tNOT_FOUND")
PY
