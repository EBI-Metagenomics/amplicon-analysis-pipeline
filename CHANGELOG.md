# EBI-Metagenomics/amplicon-analysis-pipeline: Changelog

## Unreleased

- Single-end runs now publish `primer-identification/rev_primers.fasta` in the primer's own orientation, as paired-end runs do. Previously it was reverse-complemented. Cutadapt primer files (forward, reverse and their reverse complements, padded and anchored) are now generated from it by the new `PREP_CUTADAPT_PRIMERS` step, controlled by `--bases_allowed_before_primer` (default 25).
- Cutadapt now trims primers in two passes: each read's own primer from its 5' end, then the other primer (reverse-complemented) from its 3' end, so paired-end reads that run through a short amplicon no longer keep the other primer and ASVs no longer start or end with primers. Both passes are anchored within `--bases_allowed_before_primer` of the read end. The 5' pass uses `-e 0.2` (previously `-n 2` with the default 0.1); the 3' pass uses `-e 0.1 -O 7`, so chance matches of a few primer bases at read ends are not trimmed, at the cost of leaving read-through primer fragments shorter than 7 nt. Each run now writes `*.first_pass.cutadapt.json` and `*.second_pass.cutadapt.json` instead of `*.cutadapt.json`.

## v6.2 - [2026/09/08]

- The identification of ITS-containing sequences has been extended to include the 5.8S rRNA region. Previously, sequences matching SSU and LSU Rfam models had their corresponding regions masked before the remaining sequence was queried against ITS reference databases. The same approach is now applied to 5.8S regions, if provided, to allow the annotation of datasets targeting 5.8S+ITS.
- Solved a bug that prevented the annotation of LSU sequences.
- Introduced a `target` flag for reference databases, allowing the pipeline to identify which databases to apply to each sequence type, rather than launching all databases on all sequences. The pipeline uses this field to:
  - Run SSU databases against SSU sequences only
  - Run LSU databases against LSU sequences only
  - Avoid running every reference database against every sequence type
  - ITS database searches remain an exception: they are run for all relevant sequences regardless of the configured target. This is necessary because ITS sequences are identified by masking SSU and LSU, including 5.8S regions.

## v6.1 - [2026/04/23]

- Update to publish all ASVs even if they do not have a taxonomic assignment.
- Added additional dada2 summary stats including automatically-chosen truncation points and read counts at intermediate filtering steps.
- Added flexibility for defining what MapSeq databases are used: instead of keeping a strict definition matching MGnify-specific database structure (e.g. UNITE and ITSoneDB for ITS), the new implementation allows for a different number of databases to be selected (e.g. only one db, which could be different from UNITE or ITSoneDB as long as it's similarly formatted). This also allows the usage of user-specific customised databases. 
- To allow greater flexibility when performing bespoke analyses, we enabled three different merging modes in DADA2:
  1. Merge (current behaviour)
  2. Merge-with-gaps (keep ASVs that don't merge with an overlap, instead inserting a gap)
  3. Keep separate and concatenate ASVs deriving from different paired-ends


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

