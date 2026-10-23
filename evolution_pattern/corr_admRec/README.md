# Recombination–Admixture Correlation

This project evaluates the relationship between recombination rate (cM/Mb) and *Zea mays* ssp. *mexicana* admixture dosage. Analyses use B73 RefGen_v4 coordinates and cover both population-average admixture and per-individual estimates.

## Project relationship

This directory is maintained as a standalone local Git repository during analysis. Its canonical publication location is `mexicana_introgression/evolution_pattern/corr_admRec/` in the parent GitHub repository. Do not push this standalone repository directly to the parent repository root.

## Repository layout

- `scripts/`: R analysis programs, LSF job files, and the pipeline command record.
- `data/raw/`: read-only source data; source datasets currently reside outside this repository.
- `data/processed/`: reproducible intermediate data generated from raw inputs.
- `results/<date>_<description>/`: run-specific outputs. Track only `RUN.md` and `SUMMARY.md` in Git.
- `log/`: local scheduler and program logs.
- `old_canbedeleted/`: archived legacy results; never add this directory to Git.

## Workflow

Jobs run on the NCPGR LSF cluster. Review paths and parameters in `scripts/pipeline_corr_RecAdmix.sh` and the referenced `.lsf` files before submission. The current jobs load `Anaconda3/4.3.1` and activate `R4.4_env`.

From `scripts/`, the principal stages are:

```bash
bsub < addRec_to_admMap_perChr.lsf
bsub < plot_corr_recVSadmix_allChrs.lsf
bsub < correlation_decileScatter_recVSadmix_allChrs.lsf
```

Do not submit jobs until input paths, sample lists, reference coordinates, output directory, queue, CPUs, and memory have been verified.

## Reproducibility

Each new run must use `results/<date>_<description>/` and include:

- `RUN.md`: commands, input MD5 values, software versions, parameters, runtime, job IDs, and output inventory.
- `SUMMARY.md`: sample size, statistical tests, multiple-testing correction, thresholds, effect sizes, key findings, and limitations.

Randomized analyses must use seed `20261023`. Generated data, plots, logs, and biological inputs remain local and are excluded from Git.

## Publishing to the parent repository

After the analysis is complete:

1. Commit the reviewed scripts, this README, and run summary documents locally.
2. Update a separate, clean checkout of the `mexicana_introgression` repository.
3. Copy only locally tracked files into `evolution_pattern/corr_admRec/`, preserving their relative paths.
4. Review the parent repository's status and diff before committing there.
5. Push from the parent repository only after confirming that no data, generated figures, logs, or archived results are included.

The standalone repository intentionally has no remote configured; the parent-repository checkout is the publication boundary.
