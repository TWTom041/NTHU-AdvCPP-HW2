#!/usr/bin/env bash
# Fetch the public RAPTOR code (Bronzwaer, Davelaar et al. 2018, GPLv2) at a pinned commit.
set -euo pipefail
cd "$(dirname "$0")/.."
DEST=third_party/raptor
COMMIT=08cb9a2bba526dc7f0ee91e59ff7e178d0e709a1
if [ ! -d "$DEST/.git" ]; then
  git clone https://github.com/tbronzwaer/raptor.git "$DEST"
fi
git -C "$DEST" checkout -q "$COMMIT"
echo "RAPTOR at $(git -C "$DEST" rev-parse HEAD)"
