#!/usr/bin/env bash
set -euo pipefail

# Update a consumer module to a published bamdriver-go version.
# Usage:
#   ./scripts/update_consumer.sh /path/to/consumer v0.1.0

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <consumer_dir> <version-tag>"
  exit 1
fi

CONSUMER_DIR="$1"
VERSION="$2"

if [[ ! "$VERSION" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "error: version must match vMAJOR.MINOR.PATCH"
  exit 1
fi

if [[ ! -f "$CONSUMER_DIR/go.mod" ]]; then
  echo "error: $CONSUMER_DIR/go.mod not found"
  exit 1
fi

cd "$CONSUMER_DIR"

# Update requirement
if rg -q "github.com/rainoffallingstar/bamdriver-go" go.mod; then
  go mod edit -require=github.com/rainoffallingstar/bamdriver-go@"$VERSION"
else
  go mod edit -require=github.com/rainoffallingstar/bamdriver-go@"$VERSION"
fi

# Remove local replace if present
if rg -q "replace github.com/rainoffallingstar/bamdriver-go" go.mod; then
  go mod edit -dropreplace=github.com/rainoffallingstar/bamdriver-go || true
fi

go mod tidy
go test ./...

echo "[done] updated $(pwd) to bamdriver-go $VERSION"
