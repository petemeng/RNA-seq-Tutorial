# PRJDB11848 Run Snapshot

This folder is a published snapshot of one validated run, intended for web inspection and reproducibility checks.

## Included

- `metadata/`: sample sheet and metadata files used by the run
- `logs/`: download and quantification logs
- `results/`: chapter 4-9 downstream outputs
- `tool_versions.txt`: runtime tool version record
- `FILELIST.txt`: complete list of files in this snapshot
- `CHECKSUMS.sha256`: checksums for all files (except itself)

## Quick checks

```bash
# WT DEG count
tail -n +2 artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv | wc -l
# expected: 1541

# WT time-series significant genes
tail -n +2 artifacts/prjdb11848/results/ch8/time_series_LRT_WT_sig.csv | wc -l
# expected: 2927

# verify checksums
sha256sum -c artifacts/prjdb11848/CHECKSUMS.sha256
```
