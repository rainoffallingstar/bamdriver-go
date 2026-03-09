#!/usr/bin/env bash
set -euo pipefail

# Release helper for bamdriver-go.
# Usage:
#   ./scripts/release.sh v0.1.0

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <version-tag>"
  exit 1
fi

TAG="$1"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$ROOT_DIR"

if [[ ! "$TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "error: tag must match vMAJOR.MINOR.PATCH"
  exit 1
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "error: not a git repository"
  exit 1
fi

if [[ -n "$(git status --porcelain)" ]]; then
  echo "error: working tree is not clean"
  git status --short
  exit 1
fi

go mod tidy
go test ./...

git tag "$TAG"

echo "[done] created tag $TAG"
echo "next: git push origin main --tags"
