#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COMMIT="fd22dec943a3a9cc70643f079c40a8baf13a8e62"
DEST_DIR="$ROOT/.external/NanoSeq_${COMMIT}"
DEST="$DEST_DIR/extract_tags.py"
URL="https://raw.githubusercontent.com/cancerit/NanoSeq/${COMMIT}/python/extract_tags.py"

mkdir -p "$DEST_DIR"
if ! command -v curl >/dev/null 2>&1; then
  echo "ERROR: curl is required to fetch the pinned NanoSeq file." >&2
  exit 2
fi
curl --fail --location "$URL" --output "$DEST"
chmod +x "$DEST"
cat <<MSG
Fetched pinned NanoSeq extract_tags.py:
  $DEST
Source:
  $URL
Upstream license: GNU AGPL v3 or later (retain upstream notices).

The historical scNanoSeq preprocessing script expects extract_tags.py to be
available in its execution/script directory. Copy or symlink this exact pinned
file into that environment, or update the path explicitly and record the change.
MSG
