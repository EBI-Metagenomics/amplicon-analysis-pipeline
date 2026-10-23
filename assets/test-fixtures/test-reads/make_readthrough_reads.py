#!/usr/bin/env python3
"""Generate synthetic 16S V4 (515F/806R) reads where the read runs through the amplicon.

The amplicon is 292 bp with primers (253 bp without), so 2x300 reads end in the
other primer, reverse-complemented, followed by adapter. Used to test that
primers are trimmed from both ends of each read.

Read classes, written as `class=<name>` in each read's header:
    through   insert shorter than the read: ... RC(primer) + adapter + G fill
    edge      read ends exactly at the insert edge, as fastp's overlap trimming leaves it
    partial   read ends partway into the RC primer, as quality trimming leaves it
    long      insert longer than the read (150 nt reads): no read-through
    internal  long read with an RC(806R) motif placed 80 nt in: must not be cut there

Outputs (paths relative to the repo root):
    assets/test-fixtures/test-reads/READTHROUGH_{1,2}.fastq.gz     class `through`, end-to-end test (6000 pairs)
    modules/ebi-metagenomics/cutadapt/tests/data/readthrough_{1,2}.fastq.gz   all classes, PE
    modules/ebi-metagenomics/cutadapt/tests/data/readthrough_se.fastq.gz      all classes, SE

Usage: python3 assets/test-fixtures/test-reads/make_readthrough_reads.py  (from the repo root)
"""

import gzip
import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SSU = REPO / "assets/test-fixtures/test-ref-dbs/SILVA-SSU/SILVA-SSU_small.fasta"
E2E_OUT = REPO / "assets/test-fixtures/test-reads"
MODULE_OUT = REPO / "modules/ebi-metagenomics/cutadapt/tests/data"

# Concrete primer sequences as they occur in both SILVA-SSU_small templates
# (degenerate forms: 515F GTGYCAGCMGCCGCGGTAA, 806R GGACTACNVGGGTWTCTAAT).
FWD = "GTGCCAGCAGCCGCGGTAA"
REV = "GGACTACCAGGGTATCTAAT"
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
READ_LEN = 300
LONG_READ_LEN = 150
INTERNAL_POS = 80
PARTIAL_CUTS = (5, 8, 12)  # nt of the RC primer missing from the read end
ERROR_RATE = 0.002  # end-to-end fixture only, so DADA2's learnErrors has errors to learn

COMP = str.maketrans("ACGT", "TGCA")


def rc(seq):
    return seq.translate(COMP)[::-1]


def templates():
    seqs, name = {}, None
    for line in SSU.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line.strip().upper().replace("U", "T")
    amplicons = []
    for seq in seqs.values():
        start = seq.index(FWD)
        end = seq.index(rc(REV)) + len(REV)
        amplicons.append(seq[start:end])
    # A third template: the first with 8 substitutions inside the amplicon (~3%),
    # away from the primers, so DADA2 has a close variant to resolve.
    variant = list(amplicons[0])
    for pos in range(40, 240, 25):
        variant[pos] = {"A": "G", "G": "A", "C": "T", "T": "C"}[variant[pos]]
    amplicons.append("".join(variant))
    return amplicons


def fill(seq, adapter, length):
    return (seq + adapter + "G" * length)[:length]


def mutate(seq, rng):
    return "".join(
        rng.choice([b for b in "ACGT" if b != base]) if rng.random() < ERROR_RATE else base
        for base in seq
    )


def read_pair(amplicon, cls, cut=0):
    """Return (R1, R2) for one pair of the given class."""
    r1_insert, r2_insert = amplicon, rc(amplicon)
    if cls == "through":
        return fill(r1_insert, ADAPTER_R1, READ_LEN), fill(r2_insert, ADAPTER_R2, READ_LEN)
    if cls == "edge":
        return r1_insert, r2_insert
    if cls == "partial":
        return r1_insert[:-cut], r2_insert[:-cut]
    if cls == "long":
        return r1_insert[:LONG_READ_LEN], r2_insert[:LONG_READ_LEN]
    if cls == "internal":
        r1 = r1_insert[:INTERNAL_POS] + rc(REV) + r1_insert[INTERNAL_POS + len(REV):]
        r2 = r2_insert[:INTERNAL_POS] + rc(FWD) + r2_insert[INTERNAL_POS + len(FWD):]
        return r1[:LONG_READ_LEN], r2[:LONG_READ_LEN]
    raise ValueError(cls)


def fastq_record(name, n, mate, cls, seq, rng):
    # The /1 or /2 suffix must end the line (FASTQSUFFIXHEADERCHECK), so the class goes before it
    header = f"@{name}.{n} class={cls} M00000:1:000000000-AAAAA:1:1101:{1000 + n}:{2000 + n}/{mate}"
    qual = "".join(rng.choice("FGHI") for _ in seq)  # Q37-40: trunc_len_automation won't truncate
    return f"{header}\n{seq}\n+\n{qual}\n"


def write(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 keeps the gzip bytes identical between regenerations
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
        fh.write("".join(records).encode())


def main():
    rng = random.Random(42)
    amplicons = templates()

    # End-to-end: class `through` only, with sequencing errors. 6000 pairs: classify_var_regions
    # needs >= 5000 sequences per region (MIN_SEQ_COUNT, mgnify-pipelines-toolkit)
    r1s, r2s, n = [], [], 0
    for _ in range(2000):
        for amplicon in amplicons:
            n += 1
            r1, r2 = read_pair(amplicon, "through")
            r1s.append(fastq_record("READTHROUGH", n, 1, "through", mutate(r1, rng), rng))
            r2s.append(fastq_record("READTHROUGH", n, 2, "through", mutate(r2, rng), rng))
    write(E2E_OUT / "READTHROUGH_1.fastq.gz", r1s)
    write(E2E_OUT / "READTHROUGH_2.fastq.gz", r2s)

    # Module: every class, error-free so trimmed lengths are exact
    cases = [("through", 0), ("edge", 0), ("long", 0), ("internal", 0)]
    cases += [("partial", cut) for cut in PARTIAL_CUTS]
    r1s, r2s, ses, n = [], [], [], 0
    for amplicon in amplicons:
        for cls, cut in cases:
            for _ in range(2):
                n += 1
                r1, r2 = read_pair(amplicon, cls, cut)
                r1s.append(fastq_record("readthrough", n, 1, cls, r1, rng))
                r2s.append(fastq_record("readthrough", n, 2, cls, r2, rng))
                ses.append(fastq_record("readthrough_se", n, 1, cls, r1, rng))
    write(MODULE_OUT / "readthrough_1.fastq.gz", r1s)
    write(MODULE_OUT / "readthrough_2.fastq.gz", r2s)
    write(MODULE_OUT / "readthrough_se.fastq.gz", ses)


if __name__ == "__main__":
    main()
