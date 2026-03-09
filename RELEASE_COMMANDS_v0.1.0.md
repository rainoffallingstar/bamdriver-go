# Release Commands (v0.1.0)

This runbook assumes three local repos:

- `bamdriver-go`
- `xenofilter-go`
- `Paireads`

## 0) Preconditions

```bash
# Ensure clean trees
cd /home/fallingstar10/xdxtools/bamdriver-go && git status --short
cd /home/fallingstar10/xdxtools/xenofilter-go && git status --short
cd /home/fallingstar10/xdxtools/Paireads && git status --short
```

## 1) Release bamdriver-go v0.1.0

```bash
cd /home/fallingstar10/xdxtools/bamdriver-go

# Optional: initialize standalone repo if needed
./scripts/bootstrap_repo.sh <REMOTE_URL>

# Commit release content
git add .
git commit -m "release: v0.1.0"

# Validate + tag
./scripts/release.sh v0.1.0

# Publish
git push origin main --tags
```

## 2) Upgrade xenofilter-go to v0.1.0

```bash
cd /home/fallingstar10/xdxtools/bamdriver-go
./scripts/update_consumer.sh ../xenofilter-go v0.1.0

cd /home/fallingstar10/xdxtools/xenofilter-go
git add go.mod go.sum internal/bgzip internal/bamnative cmd/test_elprep internal/bam internal/bgzip
git commit -m "refactor: use bamdriver-go v0.1.0"
git push
```

## 3) Upgrade Paireads to v0.1.0

```bash
cd /home/fallingstar10/xdxtools/bamdriver-go
./scripts/update_consumer.sh ../Paireads v0.1.0

cd /home/fallingstar10/xdxtools/Paireads
git add go.mod go.sum bamnative internal/bgzip
git commit -m "refactor: use bamdriver-go v0.1.0"
git push
```

## 4) Post-release verification

```bash
cd /home/fallingstar10/xdxtools/xenofilter-go && go test ./...
cd /home/fallingstar10/xdxtools/Paireads && go test ./...
```

If `samtools` is available:

```bash
samtools quickcheck <output.bam>
samtools view -H <output.bam> >/dev/null
samtools view -c <output.bam> >/dev/null
```

## 5) Rollback (if needed)

```bash
# In consumers, pin back to local replace (or previous tag) and retest
cd /home/fallingstar10/xdxtools/xenofilter-go
# edit go.mod to previous working state, then:
go mod tidy && go test ./...

cd /home/fallingstar10/xdxtools/Paireads
# edit go.mod to previous working state, then:
go mod tidy && go test ./...
```
