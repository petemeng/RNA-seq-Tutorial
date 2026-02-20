---
name: rnaseq-arabidopsis-doc-maintainer
description: Maintain and expand this repository RNA-seq tutorial with a single Arabidopsis thaliana Mock/flg22 storyline. Use when editing docs chapters, unifying naming, fixing code snippets, checking consistency, and preparing documentation updates for debugging or revision.
---

# RNA-seq Arabidopsis Doc Maintainer

## Fixed Scope

- Use Arabidopsis thaliana examples only.
- Use condition names `Mock` and `flg22` only.
- Prefer tutorial run examples `SRR8694017/18/19/20`.
- Keep TAIR10 references:
  - `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa`
  - `Arabidopsis_thaliana.TAIR10.58.gtf`
  - `Arabidopsis_thaliana.TAIR10.cdna.all.fa`

## Chapter Contract

- Keep `Why over How` explanation.
- Keep at least one Pitfall block per chapter.
- Keep runnable Bash or R code blocks.
- Keep a chapter-end checklist.

## Technical Defaults

- Upstream stack: `fastp + MultiQC + STAR (2-pass) + Salmon`.
- Differential analysis: `DESeq2` with default contrast `flg22 vs Mock`.
- Enrichment: `clusterProfiler + org.At.tair.db`.
- KEGG organism code: `ath`.

## Naming and Code Rules

- Keep `Mock/flg22`; avoid `Control/Treated`.
- Do not reintroduce human reference paths.
- Run `resultsNames(dds)` before extracting DESeq2 coefficients.

## Workflow

1. Scan for inconsistent terms and stale examples.
2. Edit target chapters while preserving structure.
3. Re-scan to confirm no regressions.
4. If user asks to publish, commit docs-only changes and push.

## Verification Commands

```bash
# Negative checks: should not exist
rg -n "GRCh38|GENCODE|Homo sapiens|org\.Hs\.eg\.db|\bhsa\b|Treated|Control|SRR1039508|SRR1039509" docs README.md

# Positive checks: should exist
rg -n "Arabidopsis|TAIR10|org\.At\.tair\.db|\bath\b|Mock|flg22|SRR8694017|SRR8694019|SRR8694020" docs README.md
```

## Guardrails

- Do not modify `.claude/settings.local.json`.
- Do not reorder `mkdocs.yml` nav unless explicitly requested.
- If adding chapters, sync entries in both `docs/index.md` and `mkdocs.yml`.
