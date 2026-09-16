# EBI-Metagenomics/amplicon-analysis-pipeline: Changelog

## v6.2 - [2026/09/08]

- The identification of ITS-containing sequences has been extended to include the 5.8S rRNA region. Previously, sequences matching SSU and LSU Rfam models had their corresponding regions masked before the remaining sequence was queried against ITS reference databases. The same approach is now applied to 5.8S regions, if provided, to allow the annotation of datasets targeting 5.8S+ITS.
- Solved a bug that prevented the annotation of LSU sequences.
- Introduced a target flag for reference databases, allowing the pipeline to identify which databases to apply to each sequence type, rather than launching all databases on all sequences. The pipeline uses this field to:
  - Run SSU databases against SSU sequences only
  - Run LSU databases against LSU sequences only
  - Avoid running every reference database against every sequence type
  - ITS database searches remain an exception: they are run for all relevant sequences regardless of the configured target. This is necessary because ITS sequences are identified by masking SSU, LSU, or 5.8S regions.

## v6.1 - [2026/04/23]

- Update to publish all ASVs even if they do not have a taxonomic assignment.
- Added additional dada2 summary stats including automatically-chosen truncation points and read counts at intermediate filtering steps.
- Added flexibility for defining what MapSeq databases are used and how dada2 is run.


## v6.0 - [2025/10/31]

Initial release of v6 amplicon-analysis-pipeline. Re-implements all of the existing features from v5.0:

- Reads quality control
- rRNA sequence extraction using [Infernal/cmsearch](https://github.com/EddyRivasLab/infernal/tree/master)
- Closed-reference-based taxonomic classification and visualisation of rRNA using [MAPseq](https://github.com/meringlab/MAPseq) and [Krona](https://github.com/marbl/Krona)

v6.0 also contains multiple significant changes:

- Refactoring from CWL to [Nextflow](https://www.nextflow.io/) for pipeline definition
- Simplification of reads quality control using [fastp](https://github.com/OpenGene/fastp)
- Automatic amplified region inference for 16S and 18S rRNA
- Automatic primer identification, trimming, and validation
- Addition of Amplicon Sequence Variant (ASV) calling using [DADA2](https://benjjneb.github.io/dada2/index.html)
- Taxonomic classification and visualisation of ASVs using [MAPseq](https://github.com/meringlab/MAPseq) and [Krona](https://github.com/marbl/Krona) to complement the existing closed-reference analysis
- Addition of [PR2](https://pr2-database.org/) as a reference database
- Updating of existing reference databases ([SILVA](https://www.arb-silva.de/), [UNITE](https://unite.ut.ee/), [ITSoneDB](https://itsonedb.cloud.ba.infn.it), [Rfam](https://rfam.org/))

