# bamdriver-go

Shared pure-Go BAM/BGZF driver package extracted from `xenofilter-go` and `Paireads`.

## Packages

- `pkg/bgzip`: BGZF reader/writer and virtual offset helpers.
- `pkg/bamnative`: BAM reader/writer/sort/index, plus FASTA/NM helpers used by xenofilter.

## Local Development

In consumers, use:

- `require github.com/rainoffallingstar/bamdriver-go v0.0.0`
- `replace github.com/rainoffallingstar/bamdriver-go => ../bamdriver-go`

For release, replace `v0.0.0` with a tagged version and remove `replace`.

## Scripts

- `scripts/bootstrap_repo.sh [remote_url]`: initialize standalone git repo and optional origin.
- `scripts/release.sh vX.Y.Z`: run tidy/tests and create a release tag.
- `scripts/update_consumer.sh /path/to/consumer vX.Y.Z`: switch a consumer from local `replace` to a published tag.
