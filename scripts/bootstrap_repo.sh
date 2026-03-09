#!/usr/bin/env bash
set -euo pipefail

# Initialize bamdriver-go as a standalone git repo and optionally set remote.
# Usage:
#   ./scripts/bootstrap_repo.sh [remote_url]

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOTE_URL="${1:-}"

cd "$ROOT_DIR"

if [[ -d .git ]]; then
  echo "[info] .git already exists in $ROOT_DIR; skipping init"
else
  git init
  git add .
  git commit -m "chore: bootstrap bamdriver-go module"
fi

if [[ -n "$REMOTE_URL" ]]; then
  if git remote get-url origin >/dev/null 2>&1; then
    git remote set-url origin "$REMOTE_URL"
  else
    git remote add origin "$REMOTE_URL"
  fi
  echo "[info] origin -> $REMOTE_URL"
fi

echo "[done] bootstrap complete"
