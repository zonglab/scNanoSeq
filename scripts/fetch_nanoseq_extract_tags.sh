#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DEST="$ROOT/2_python_script/extract_tags.py"
URL="https://raw.githubusercontent.com/cancerit/NanoSeq/fd22dec943a3a9cc70643f079c40a8baf13a8e62/python/extract_tags.py"

printf 'Fetching pinned NanoSeq extract_tags.py\n  source: %s\n  dest:   %s\n' "$URL" "$DEST"

if command -v curl >/dev/null 2>&1; then
    curl -L --fail --silent --show-error "$URL" -o "$DEST"
elif command -v wget >/dev/null 2>&1; then
    wget -q "$URL" -O "$DEST"
else
    echo "ERROR: curl or wget is required." >&2
    exit 1
fi

chmod +x "$DEST"
python3 "$DEST" --help >/dev/null

echo "Installed pinned NanoSeq helper successfully."
echo "Keep the upstream licence/copyright header intact if you redistribute this file."
