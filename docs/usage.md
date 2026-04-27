# ebi-metagenomics/amplicon-analysis-pipeline: Usage

> [!NOTE]
> This pipeline follows some of the conventions and uses tools from the nf-core community, but it is not part of the nf-core project.

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The `ebi-metagenomics/amplicon-analysis-pipeline` analyses amplicon sequencing reads from single-end, paired-end, or interleaved FASTQ input. It performs reads QC, amplified-region inference, primer handling, ASV calling, and taxonomic classification.

## Samplesheet Input

You need to create a samplesheet before running the pipeline and pass it with `--input`. The file must be a comma-separated file with a header row.

```bash
--input ./samplesheet.csv
```

An example samplesheet is provided in [assets/samplesheet.csv](../assets/samplesheet.csv).

### Supported layouts

The samplesheet must contain these columns:

| Column       | Description                                                                                            |
| ------------ | ------------------------------------------------------------------------------------------------------ |
| `sample`     | Custom sample name. Must not contain spaces.                                                           |
| `fastq_1`    | Path or URL to read 1 FASTQ file. Files must end in `.fq.gz` or `.fastq.gz`.                           |
| `fastq_2`    | Path or URL to read 2 FASTQ file for paired-end data. Leave empty for single-end or interleaved input. |
| `single_end` | Boolean flag indicating whether the input is single-end (`true`) or paired-end/interleaved (`false`).  |

Example:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,single_end
SAMPLE_PAIRED_END,/path/to/reads/SAMPLE_PAIRED_END_1.fastq.gz,/path/to/reads/SAMPLE_PAIRED_END_2.fastq.gz,false
SAMPLE_SINGLE_END,/path/to/reads/SAMPLE_SINGLE_END.fastq.gz,,true
```

### Input rules

- `sample` values are used as output labels and cannot contain spaces.
- `fastq_1` is always required.
- `fastq_2` must be empty for single-end input.
- Interleaved paired-end reads can be supplied as a single FASTQ with `single_end` set to `false`.

## Running The Pipeline

The typical command for running the pipeline is:

```bash
nextflow run ebi-metagenomics/amplicon-analysis-pipeline \
    -r main \
    -profile docker \
    --input ./samplesheet.csv \
    --outdir ./results
```

This will create the results directory specified by `--outdir`, plus the usual Nextflow work directory and log files in your working directory.

If you want to reuse the same parameters for multiple runs, use a params file instead of writing long command lines repeatedly.

```bash
nextflow run ebi-metagenomics/amplicon-analysis-pipeline \
    -profile docker \
    -params-file params.yaml
```

### Reproducibility

It is a good idea to pin the pipeline version when running analysis on important data. Use `-r` with a release tag so Nextflow runs a fixed version of the pipeline code.

```bash
nextflow run ebi-metagenomics/amplicon-analysis-pipeline -r <release-tag> -profile docker --input ./samplesheet.csv --outdir ./results
```

You can refresh the locally cached pipeline code with:

```bash
nextflow pull ebi-metagenomics/amplicon-analysis-pipeline
```

### Core Nextflow Arguments

> [!NOTE]
> These options are part of Nextflow and use a single hyphen. Pipeline parameters use `--`.

#### `-profile`

Use `-profile` to choose a configuration profile for your environment. Common examples include `docker`, `singularity`, `podman`, and `apptainer`. Profiles can be combined, for example `-profile test,docker`.

If `-profile` is omitted, Nextflow will run locally and expect dependencies on the `PATH`, which is not recommended.

#### `-resume`

Use `-resume` to continue a previous run from cached results when the inputs have not changed.

#### `-c`

Use `-c` to load a Nextflow config file when you need to tune execution settings such as resource requests or module arguments.

### Pipeline-Specific Configuration

Some pipeline behavior is controlled through configuration rather than the command line:

- `params.mapseq_databases`
  - Database definitions used for taxonomy assignment.
  - Each entry should define `fasta`, `tax`, `otu`, `mscluster`, `run_otu`, `run_asv`, and `label`.
  - If `run_asv` is `true`, the entry also needs `asv_label`.
- `params.rrnas_rfam_covariance_model`
  - Path to the Rfam covariance model directory used for rRNA detection.
- `params.rrnas_rfam_claninfo`
  - Path to the Rfam clan info file used for rRNA detection.
- `params.std_primer_library`
  - Optional path to a directory containing a standard primer library for primer inference and validation.

If you are customising these settings, prefer a config file or a params file so the run is reproducible.
