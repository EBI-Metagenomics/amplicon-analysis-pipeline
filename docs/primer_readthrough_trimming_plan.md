# Trimming read-through primers at 3′ read ends — implementation plan

Oct 20, 2026 · @Tim

Builds on `~/Downloads/primer-adapter-trimming-pipelines.md` (Oct 8, 2026). Reviewed against
`da2433d` on `main`, and against envIdent `05af595`, which includes PR #12 ("Restructuring primer trimming").

## Problem

In paired-end runs, when the amplicon is shorter than the read:

- R1 keeps the reverse primer, reverse-complemented, at its 3′ end.
- R2 keeps the forward primer, reverse-complemented, at its 3′ end.

How it happens:

- fastp (`READS_QC`) overlap-trims each read back to the outer edge of the amplicon. The primers are still on the reads at that point, so the cut lands *outside* the other primer.
- `CUTADAPT` then gets `-g file:fwd -G file:rev` and `ext.args = "-n 2"` (`conf/modules.config:244`), which removes only the primer at the start of each read.
- `bin/dada2.R:209` calls `mergePairs` without `trimOverhang`, so overhangs stay in the merged sequence. ASVs can therefore start with the forward primer and end with the reverse-complemented reverse primer. This happens when `trunc_len_automation` doesn't truncate the reads first.

Single-end runs are already handled. `PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS` passes `--se`, which writes `rev_primers.fasta` reverse-complemented, and `CUTADAPT` uses it with `-a`.

## Strategy (as in envIdent PR #12)

- Write four primer files per run: `F`, `R`, `F_RC`, `R_RC`. Pad them with `N`s and anchor them with `X`. `bases_allowed_before_primer` sets the padding and defaults to 25.
- Run cutadapt twice, with `-e 0.2`:

| Pass | R1 | R2 | Removes |
| --- | --- | --- | --- |
| 1 | `-g XN{25}FWD` | `-G XN{25}REV` | The read's own primer, within 25 nt of the start |
| 2 | `-a RC(REV)N{25}X` | `-A RC(FWD)N{25}X` | The other primer, within 25 nt of the end (partial occurrences at the end included) |

Single-end runs use `-g` in pass 1 and `-a` in pass 2. `X` stops primer-like sequence in the middle of a read from truncating it. Two passes mean `-n 2` is no longer needed.

---

## Phase 0 — Test that reproduces the bug (do first, commit red)

Goal: a small, deterministic fixture whose tests **fail on current `main`** and pass after the fix.

### 0.1 Fixture generator

`assets/test-fixtures/test-reads/make_readthrough_reads.py`: stdlib only, fixed seed, committed next to its output so the fixture can be regenerated.

- **Templates:** 3–5 real 16S sequences, so that `DETECT_RNA` (cmsearch against `ribo.cm`) and MAPseq classify them. Start with the two in `test-ref-dbs/SILVA-SSU/SILVA-SSU_small.fasta`, plus a few full-length references if needed.
- **Primers:** V4, 515F `GTGYCAGCMGCCGCGGTAA` and 806R `GGACTACNVGGGTWTCTAAT`. Pick one concrete base at each degenerate position for the reads. The amplicon is about 253 bp between the primers, or about 292 bp with them.
- **Reads:** 2×300 with flat high quality (for example `I`), so that `trunc_len_automation` returns 0 and DADA2 doesn't truncate the primer away. Apply about 0.2% random substitutions so `learnErrors` has something to learn from.
- **Adapters after the RC primer:** TruSeq, `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` on R1 and `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` on R2, then `G`s to fill the read.

Pair classes, labelled in the read header so the tests can select them:

| Class | What the read looks like | Expected after the fix |
| --- | --- | --- |
| A `through` | Insert shorter than the read: … RC(primer) + adapter | Both primers and the adapter tail removed |
| B `edge` | Read ends exactly at the insert edge, as fastp leaves it | Both primers removed |
| C `partial` | Read ends partway into the RC primer (8–15 nt), as quality trimming leaves it | Partial RC primer removed |
| D `long` | Insert longer than the read: no read-through | Only the 5′ primer removed. 3′ end untouched |
| E `internal` | RC(806R)-like motif placed about 100 nt into the read | Read **not** truncated at the motif |

Outputs:

- `READTHROUGH_{1,2}.fastq.gz`: class A only, 2,000 pairs per template (6,000 pairs), for the end-to-end test. `classify_var_regions` needs at least 5,000 sequences per region (`MIN_SEQ_COUNT` in mgnify-pipelines-toolkit); with fewer, the run ends as `no_asvs` before any primer step.
- `modules/ebi-metagenomics/cutadapt/tests/data/readthrough_{1,2}.fastq.gz` and `readthrough_se.fastq.gz`: all classes, about 20 pairs each, for the module test.
- `modules/ebi-metagenomics/cutadapt/tests/data/{fwd,rev}_primers.fasta`: 515F and 806R, as `PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS` writes them for a paired-end run.

### 0.2 Module test: `modules/ebi-metagenomics/cutadapt/tests/main.nf.test` (new; there are no cutadapt tests in this repo yet)

- PE: input `[ [id:'readthrough', single_end:false], reads, primers ]`. Assert:
  - No R1 matches `ATTAGA[AT]ACCC[CGT].GTAGTCC` (RC 806R) within its last 30 nt.
  - No R2 matches `TTACCGCGGC[GT]GCTG[AG]CAC` (RC 515F) within its last 30 nt.
  - No read starts with either primer.
  - Classes A–C come out at amplicon length (about 253, ±2).
  - Class D's 3′ ends are unchanged.
  - Class E keeps its full length.
- SE: the same, on the single-end file. **Result on `main`: fails on class E only.** The unanchored `-a` cuts the internal reads at the RC(806R) motif (61 nt instead of 131). The `N{25}X` anchoring in Phase 2 fixes this too.
- No primers (`[]`): placeholder outputs exist and are empty.

### 0.3 End-to-end test: new test in `tests/main.nf.test`, tag `"readthrough"`

- Samplesheet: `assets/test-fixtures/samplesheets/test_samplesheet_readthrough.csv` with `READTHROUGH` (PE). Use a local path or the GitHub raw URL, following the existing sheets.
- First check that auto-detection finds 515F/806R in the fixture by looking at `READTHROUGH/primer-identification/*primers.fasta`. The test profile sets `std_primer_library = ""`. If auto-detection misses them, give this test a small standard primer library containing 515F/806R.
- Assert on the ASV FASTA under `READTHROUGH/asv/…`:
  - No ASV starts with 515F.
  - No ASV ends with RC(806R). Use the IUPAC regexes above.
  - ASV lengths are within 250–256.
- **Result on `main`: fails.** Auto-detection finds `V4T9F`/`V4R19R`, and DADA2 produces 3 ASVs, all 292 bp, each starting with 515F and ending with RC(806R). 5′ trimming happens, but the `mergePairs` overhang from the other mate puts each primer back. Committed red on `test/primer-readthrough`.

---

## Phase 1 — Primer files

1. **Drop `--se`** in `modules/local/primer_validation_classify_var_regions/main.nf`, so both layouts get primers in their own orientation. The existing snapshot shows `--se` changes only `rev_primers.fasta`: `primer_validation.tsv` has the same md5 either way. Update that module's snapshot: the SE test's `rev_primers.fasta` md5 changes. The published SE `rev_primers.fasta` changes orientation, so note it in `CHANGELOG.md`.
2. **Add `PREP_CUTADAPT_PRIMERS`**, ported from envIdent (`modules/local/prep_cutadapt_primers/main.nf`, `bin/revcomp_primers.py`, `bin/prep_cutadapt_primers.py`):
   - Input `[meta, fwd_primers.fasta, rev_primers.fasta]`. Output `[meta, F, R, F_RC, R_RC]`.
   - Use `--replace-i` (inosine → N), and `--n ${params.bases_allowed_before_primer}`.
   - Use the `mgnify-pipelines-toolkit` container, which already has Biopython, rather than adding a new one.
3. **`params.bases_allowed_before_primer = 25`**, in `nextflow.config` and `nextflow_schema.json`.
4. **Wire it into `subworkflows/local/concat_primer_cutadapt.nf`**, between `PRIMER_VALIDATION.out.validated_primers` and `cutadapt_input`. Runs without primers (`runs_without_primers`, `[meta, []]`) skip the prep step and mix back in afterwards, as they do now.
   - Check that `validated_primers` always emits both files (`fwd_` before `rev_` alphabetically). The current module indexes `primers[0]` and `primers[1]` and would break if one were missing.

## Phase 2 — Two-pass cutadapt

1. Patch `modules/ebi-metagenomics/cutadapt/main.nf` to match envIdent's version:
   - Four-file convention `[fwd, rev, fwd_rc, rev_rc]`.
   - Pass 1 writes `*.first_pass.fastq.gz`. Pass 2 writes `*.trim.fastq.gz`, so `dada2_input_preparation_function` needs no change.
   - JSON per pass: `*.first_pass.cutadapt.json` and `*.second_pass.cutadapt.json`.
   - Keep the "no primers" placeholder branch.
   - Regenerate `cutadapt.diff` with `nf-core modules patch` (ebi-metagenomics remote).
   - Leave out envIdent's `--info-file` / `primer_trim_stats.sh` summaries unless someone needs them.
2. Keep the cutadapt 4.6 container. `XN{n}…` and `…N{n}X` are supported, and keeping it keeps the diff small. Bumping to 5.2 to match envIdent can be a separate change.
3. In `conf/modules.config`, change `withName: CUTADAPT` from `ext.args = "-n 2"` to `ext.args = "-e 0.2"`, matching envIdent. The default is 0.1. Note in the PR that both passes now use 0.2.
4. Behaviour changes to call out in the PR:
   - The 5′ pass is now anchored within 25 nt (`XN{25}`). Reads whose primer starts further in are no longer trimmed, and primer-like hits deeper in a read no longer cut off everything before them.
   - Reads get shorter by one primer length. `filterAndTrim` drops reads shorter than `truncLen`, so check read retention (Phase 4).

Optional: upstream the two-pass module to `ebi-metagenomics` nf-modules so AAP and envIdent can drop their patches.

## Phase 3 — Downstream of the JSON change

1. **MultiQC input** (`workflows/pipeline.nf:524–545`): `cutadapt` becomes a list of two JSONs per run. Flatten `final_inputs` and keep the `if (!cutadapt)` guard working with an empty list.
2. **`assets/multiqc_config.yml`**: copy envIdent's `cutadapt_first_pass` / `cutadapt_second_pass` module entries, with `fn` patterns `*first_pass.cutadapt.json` and `*second_pass.cutadapt.json`. Otherwise MultiQC sees two cutadapt reports with the same sample name and keeps only one.
3. **Publishing**: the `CUTADAPT` publish pattern `*.json` already picks up both files.
4. **`docs/output.md`**: document the two JSONs.

## Phase 4 — Verification

1. **Phase 0 tests pass.** Module test classes A–E, the SE regression guard, and the end-to-end `readthrough` ASV assertions.
2. **Existing snapshots** (`tests/main.nf.test`, tag `mixed_amp`):
   - Expect changes only to cutadapt JSONs (now two per run), the SE `rev_primers.fasta`, MultiQC, and `MOCKDOUBLEAMP` outputs (about 99 of its 12,000 pairs show read-through).
   - `ERR4334351` is a V3–V4 amplicon (about 465 bp) read at 2×301, so it never reads through. Its ASVs and count tables should be **byte-identical**. A change there points to the 5′ anchoring or `-e 0.2`; explain it before updating the snapshot.
3. **Real runs.** Run 3–5 MGnify V4 runs (2×250 and 2×300) on `main` and on the branch. Compare:
   - the fraction of ASVs ending in RC(reverse primer) or starting with the forward primer;
   - the number of ASVs (expect some to collapse into one);
   - reads kept at `filterAndTrim` and `mergePairs`, from the DADA2 report.
4. **Past results (sizing).** Scan a sample of published AAP ASV FASTAs for primers at ASV ends, using the same regexes per primer pair. This shows how much past output is affected.

## Out of scope / noted

- Moving fastp after cutadapt. It would catch most of this without new primer files, but it means splitting `READS_QC` and gives less visible reporting. See the Oct 8 doc.
- `trimOverhang=TRUE` in `mergePairs`. It would hide what's left of the problem in ASVs, but it doesn't fix the reads. Add it only if Phase 4 still finds primer ends.
- `dada2_input_preparation_function` treats `reads[0]` as fastp and `reads[1]` as cutadapt after `.mix().groupTuple()`. That order isn't guaranteed, so it should get its own issue.
