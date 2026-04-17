# MGnify amplicon analysis pipeline

This repository contains the v6.0 [MGnify](https://www.ebi.ac.uk/metagenomics) amplicon analysis pipeline. It is, first and foremost, a refactor of the existing [v5.0 amplicon analysis pipeline](https://github.com/EBI-Metagenomics/pipeline-v5), replacing CWL with [Nextflow](https://www.nextflow.io/) as its workflow management system. This pipeline re-implements all [existing closed-reference v5.0 features](https://docs.mgnify.org/src/docs/analysis.html#amplicon-analysis-pipeline), and makes multiple significant changes and additions.

![V6 Schema](assets/v6_amplicon_schema.png)

## Pipeline description

### Features

The amplicon analysis pipeline v6.0 re-implements all of the existing features from v5.0:

- Reads quality control
- rRNA sequence extraction using [Infernal/cmsearch](https://github.com/EddyRivasLab/infernal/tree/master)
- Closed-reference-based taxonomic classification and visualiation of rRNA using [MAPseq](https://github.com/meringlab/MAPseq) and [Krona](https://github.com/marbl/Krona)

The amplicon analysis pipeline v6.0 also contains multiple significant changes:

- Refactoring from CWL to [Nextflow](https://www.nextflow.io/) for pipeline definition
- Simplification of reads quality control using [fastp](https://github.com/OpenGene/fastp)
- Automatic amplified region inference for 16S and 18S rRNA
- Automatic primer identification, trimming, and validation
- Addition of Amplicon Sequence Variant (ASV) calling using [DADA2](https://benjjneb.github.io/dada2/index.html)
- Taxonomic classification and visualisation of ASVs using [MAPseq](https://github.com/meringlab/MAPseq) and [Krona](https://github.com/marbl/Krona) to complement the existing closed-reference analysis
- Addition of [PR2](https://pr2-database.org/) as a reference database
- Updating of existing reference databases ([SILVA](https://www.arb-silva.de/), [UNITE](https://unite.ut.ee/), [ITSoneDB](https://itsonedb.cloud.ba.infn.it), [Rfam](https://rfam.org/))

### Valid amplicons

At this stage, the only sequence amplicons that this pipeline is built for are:

| Amplicon | Closed-reference analysis | ASV analysis |
| :------: | :-----------------------: | :----------: |
|   16S    |             ✓             |      ✓       |
|   18S    |             ✓             |      ✓       |
|   LSU    |             ✓             |      ✗       |
|   ITS    |             ✓             |      ✗       |

### Tools

| Tool                                                                                            | Version | Purpose                                                |
| ----------------------------------------------------------------------------------------------- | ------- | ------------------------------------------------------ |
| [fastp](https://github.com/OpenGene/fastp)                                                      | 1.0.1   | Read quality control                                   |
| [SeqFu](https://github.com/telatin/seqfu2)                                                      | 1.20.3  | FASTQ sanity checking                                  |
| [seqtk](https://github.com/lh3/seqtk)                                                           | 1.4     | FASTQ file manipulation                                |
| [SeqKit](https://bioinf.shenwei.me/seqkit/)                                                     | 2.9.0   | FASTQ file manipulation                                |
| [easel](https://github.com/EddyRivasLab/easel)                                                  | 0.49    | FASTA file manipulation                                |
| [bedtools](https://bedtools.readthedocs.io/en/latest/)                                          | 2.30.0  | FASTA sequence masking                                 |
| [Infernal/cmsearch](https://github.com/EddyRivasLab/infernal/tree/master)                       | 1.1.5   | rRNA sequence searching                                |
| [cmsearch_tblout_deoverlap](https://github.com/nawrockie/cmsearch_tblout_deoverlap/tree/master) | 0.09    | Deoverlapping of cmsearch results                      |
| [MAPseq](https://github.com/meringlab/MAPseq)                                                   | 2.1.1b  | Reference-based taxonomic classification of rRNA       |
| [Krona](https://github.com/marbl/Krona)                                                         | 2.8.1   | Krona chart visualisation                              |
| [cutadapt](https://cutadapt.readthedocs.io/en/stable/)                                          | 4.6     | Primer trimming                                        |
| [R](https://www.r-project.org/)                                                                 | 4.3.3   | R programming language (runs DADA2)                    |
| [DADA2](https://benjjneb.github.io/dada2/index.html)                                            | 1.30.0  | ASV calling                                            |
| [MultiQC](https://github.com/MultiQC/MultiQC)                                                   | 1.24.1  | Result aggregation into HTML reports                   |
| [mgnify-pipelines-toolkit](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit)        | 0.1.8   | Toolkit containing various in-house processing scripts |
| [PIMENTO](https://github.com/EBI-Metagenomics/PIMENTO)                                          | 1.0.2   | Primer inference toolkit used in the pipeline          |

### Reference databases

This pipeline uses five different reference databases. The files the pipeline uses are processed from the raw files available on each database's website, for use by MAPseq and cmsearch. We provide ready-made versions of these processed files on our FTP, which you can find here:

| Reference database                            | Version | Purpose                               | Processed file paths                                                                                                                                          |
| --------------------------------------------- | ------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [SILVA](https://www.arb-silva.de/)            | 138.1   | 16S+18S+LSU rRNA database             | https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-ssu/ https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-lsu/ |
| [PR2](https://pr2-database.org/)              | 5.0     | Protist-focused 18S+16S rRNA database | https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/pr2/                                                                                      |
| [UNITE](https://unite.ut.ee/)                 | 9.0     | ITS database                          | https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/unite/                                                                                    |
| [ITSoneDB](https://itsonedb.cloud.ba.infn.it) | 1.141   | ITS database                          | https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/itsonedb/                                                                                 |
| [Rfam](https://rfam.org/)                     | 14.10   | rRNA covariance models                | https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/rfam/                                                                                     |

> [!NOTE]  
> The preprocessed databases are generated with the [Microbiome Informatics reference-databases-preprocessing-pipeline](https://github.com/EBI-Metagenomics/reference-databases-preprocessing-pipeline).

## How to run

### Requirements

At the moment the only prerequisites for running it are Nextflow and [Docker](https://www.docker.com/)/[Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html), since all of the Nextflow processes use pre-built containers.

### Input shape

The input data for the pipeline is amplicon sequencing reads (either paired-end or single-end) in the form of FASTQ files. These files should be specified using a `.csv` samplesheet file with this format:

```
sample,fastq_1,fastq_2,single_end
SRR9674618,/path/to/reads/SRR9674618.fastq.gz,,true
SRR17062740,/path/to/reads/SRR17062740_1.fastq.gz,/path/to/reads/SRR17062740_2.fastq.gz,false
```

### Defining reference MAPseq databases

MAPseq databases for assigning taxonomy to sequences (i.e. reads or ASVs) are defined in a Nextflow config file as parameters. These definitions must adhere to a structure which describes

1. filepaths
2. flags determining how the database will be used by the pipeline
3. labels for naming results files

The structure is as follows:

- `params.mapseq_databases`
  - `<database_name>`
    - `fasta`: filepath of MAPseq FASTA file containing reference sequences
    - `tax`: filepath of MAPseq taxonomy of reference sequences
    - `otu`: filepath of MAPseq pre-clustered OTU reference
    - `mscluster`: filepath of MAPseq pre-computed k-mer clustering reference
    - `run_otu`: true/false for whether this database is used to run OTU analysis
    - `run_asv`: true/false for whether this database is used to run ASV analysis
    - `label`: label to be used for naming OTU analysis outputs
    - `asv_label`: label to be used for naming ASV analysis outputs (Required only if `run_asv = true`)

Examples can be found in `conf/test_dbs.config`.

### Execution

You can run the current version of the pipeline on SLURM like this:

```bash
nextflow run ebi-metagenomics/amplicon-analysis-pipeline \
    -r main \
    -profile codon_slurm \
    --input /path/to/samplesheet.csv \
    --outdir /path/to/outputdir
```

If you want to run the pipeline on deeply-sequenced reads, DADA2 can become a serious bottleneck. To counter this, we embed an extra step that runs after DADA2 has already failed twice with resources that are usually enough on average samples. This will massively boost time and cpu resources DADA2 asks for.

## Outputs

### Output directory structure

Example output directory structure for one run (`ERR4334351`):

```
├── ERR4334351
│   ├── amplified-region-inference
│   │   ├── ERR4334351.16S.V3-V4.txt
│   │   └── ERR4334351.tsv
│   ├── asv
│   │   ├── 16S-V3-V4
│   │   │   └── ERR4334351_16S-V3-V4_asv_read_counts.tsv
│   │   ├── ERR4334351_asv_seqs.fasta
│   │   ├── ERR4334351_DADA2-PR2_asv_tax.tsv
│   │   ├── ERR4334351_DADA2-SILVA_asv_tax.tsv
│   │   └── ERR4334351_dada2_stats.tsv
│   ├── primer-identification
│   │   ├── ERR4334351.cutadapt.json
│   │   ├── fwd_primers.fasta
│   │   ├── rev_primers.fasta
│   │   └── ERR4334351_primer_validation.tsv
│   ├── qc
│   │   ├── ERR4334351.fastp.json
│   │   ├── ERR4334351.merged.fastq.gz
│   │   ├── ERR4334351_dada2_errors.txt
│   │   ├── ERR4334351_multiqc_report.html
│   │   ├── ERR4334351_seqfu.tsv
│   │   └── ERR4334351_suffix_header_err.json
│   ├── sequence-categorisation
│   │   ├── ERR4334351_SSU.fasta
│   │   ├── ERR4334351_SSU_rRNA_archaea.RF01959.fa
│   │   ├── ERR4334351_SSU_rRNA_bacteria.RF00177.fa
│   │   └── ERR4334351.tblout.deoverlapped
│   └── taxonomy-summary
│       ├── DADA2-PR2
│       │   ├── ERR4334351_16S-V3-V4_DADA2-PR2_asv_krona_counts.txt
│       │   ├── ERR4334351_16S-V3-V4.html
│       │   └── ERR4334351_DADA2-PR2.mseq
│       ├── DADA2-SILVA
│       │   ├── ERR4334351_16S-V3-V4_DADA2-SILVA_asv_krona_counts.txt
│       │   ├── ERR4334351_16S-V3-V4.html
│       │   └── ERR4334351_DADA2-SILVA.mseq
│       ├── PR2
│       │   ├── ERR4334351.html
│       │   ├── ERR4334351_PR2.mseq
│       │   ├── ERR4334351_PR2.tsv
│       │   └── ERR4334351_PR2.txt
│       └── SILVA-SSU
│           ├── ERR4334351.html
│           ├── ERR4334351_SILVA-SSU.mseq
│           ├── ERR4334351_SILVA-SSU.tsv
│           └── ERR4334351_SILVA-SSU.txt
├── pipeline_info
│   ├── execution_report_2025-03-25_14-13-55.html
│   ├── execution_timeline_2025-03-25_14-13-55.html
│   ├── execution_trace_2025-03-25_14-13-55.txt
│   ├── pipeline_dag_2025-03-25_14-13-55.html
│   └── software_versions.yml
├── bco.json
├── study_multiqc_report.html
├── qc_passed_runs.csv
├── qc_failed_runs.csv
└── primer_validation_summary.json
```

For a more detailed description of the different output files, see the [OUTPUTS_DESCRIPTION.md](https://github.com/EBI-Metagenomics/amplicon-pipeline/blob/main/OUTPUTS_DESCRIPTION.md) file.

## Citations

This pipeline uses code developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
