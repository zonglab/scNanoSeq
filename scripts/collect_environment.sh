#!/usr/bin/env bash
set -u

printf 'scNanoSeq environment capture\n'
printf 'timestamp_utc: '; date -u +'%Y-%m-%dT%H:%M:%SZ'
printf 'hostname: '; hostname 2>/dev/null || true
printf 'uname: '; uname -a 2>/dev/null || true
if command -v lsb_release >/dev/null 2>&1; then lsb_release -a 2>/dev/null || true; fi
printf '\n'

version_cmd() {
  local name="$1"; shift
  printf '### %s\n' "$name"
  if command -v "$1" >/dev/null 2>&1; then
    ("$@" 2>&1 | head -n 8) || true
  else
    printf 'NOT FOUND ON PATH\n'
  fi
  printf '\n'
}

version_cmd bash bash --version
version_cmd python3 python3 --version
version_cmd bwa bwa 2>&1
version_cmd samtools samtools --version
version_cmd bcftools bcftools --version
version_cmd bedtools bedtools --version
version_cmd seqtk seqtk 2>&1
version_cmd cutadapt cutadapt --version
version_cmd picard picard --version
version_cmd java java -version
version_cmd cellphy.sh cellphy.sh --help
version_cmd R R --version

printf '### Python packages\n'
python3 - <<'PY' 2>/dev/null || true
mods = ['pysam', 'numpy', 'pandas']
for m in mods:
    try:
        mod = __import__(m)
        print(f'{m}: {getattr(mod, "__version__", "UNKNOWN")}')
    except Exception as e:
        print(f'{m}: NOT AVAILABLE ({e.__class__.__name__})')
PY

printf '\n### Hardware summary\n'
(command -v nproc >/dev/null 2>&1 && printf 'logical_cpus: %s\n' "$(nproc)") || true
(command -v free >/dev/null 2>&1 && free -h) || true
