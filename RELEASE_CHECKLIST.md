# bamdriver-go Release Checklist

## 1) Pre-release validation

Run in `bamdriver-go`:

```bash
go mod tidy
go test ./...
```

Run consumer smoke tests:

```bash
cd ../xenofilter-go && go test ./...
cd ../Paireads && go test ./...
```

## 2) Tag and publish

In `bamdriver-go` repo:

```bash
./scripts/release.sh v0.1.0
git push origin main --tags
```

## 3) Upgrade consumers

In `xenofilter-go/go.mod` and `Paireads/go.mod`:

1. Change:

```go
require github.com/rainoffallingstar/bamdriver-go v0.0.0
```

to:

```go
require github.com/rainoffallingstar/bamdriver-go v0.1.0
```

2. Remove local replace:

```go
replace github.com/rainoffallingstar/bamdriver-go => ../bamdriver-go
```

3. Run:

```bash
go mod tidy
go test ./...
```

Or from `bamdriver-go`:

```bash
./scripts/update_consumer.sh ../xenofilter-go v0.1.0
./scripts/update_consumer.sh ../Paireads v0.1.0
```

## 4) CI gate (recommended)

Add checks in each consumer pipeline:

```bash
go test ./...
# if available
samtools quickcheck <output.bam>
samtools view -H <output.bam> >/dev/null
samtools view -c <output.bam> >/dev/null
```

## 5) Post-release cleanup (optional)

After one stable release cycle:

- Remove shim/stub files that only forward to `bamdriver-go`.
- Keep package paths stable only if external users depend on them.
