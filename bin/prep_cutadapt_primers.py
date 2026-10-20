#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

"""
Add a literal string to every primer sequence in a FASTA file,
leaving headers unchanged. Optionally replace I's with N's in sequences.

If --position start: prepends "XN{N}" to the start of the sequence.
If --position end:   appends "N{N}X" to the end of the sequence.

Usage:
    python prepend_fasta.py input.fasta output.fasta
    python prepend_fasta.py input.fasta output.fasta --n 25
    python prepend_fasta.py input.fasta output.fasta --position end
    python prepend_fasta.py input.fasta output.fasta --replace-i
"""

import argparse
import sys


def modify_fasta(in_path, out_path, addition, position, replace_i):
    n_primers = 0
    with open(in_path, "r") as fin, open(out_path, "w") as fout:
        seq_chunks = []
        header = None

        def flush():
            nonlocal seq_chunks, header, n_primers
            if header is None:
                return
            seq = "".join(seq_chunks)
            if replace_i:
                seq = seq.replace("I", "N").replace("i", "n")
            if position == "start":
                seq = addition + seq
            else:
                seq = seq + addition
            fout.write(header + "\n")
            fout.write(seq + "\n")
            n_primers += 1

        for line in fin:
            line = line.strip()
            if line.startswith(">"):
                flush()
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
        flush()  # last record

    return n_primers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output FASTA file")
    parser.add_argument(
        "--n",
        type=int,
        default=10,
        help="Integer to use inside the literal string (default: 10)",
    )
    parser.add_argument(
        "--position",
        choices=["start", "end"],
        default="start",
        help="Where to add the string: 'start' prepends 'XN{N}', 'end' appends 'N{N}X' (default: start)",
    )
    parser.add_argument(
        "--replace-i",
        action="store_true",
        help="Replace any I's (and i's) in the sequences with N's (and n's)",
    )
    args = parser.parse_args()

    if args.position == "start":
        addition = f"XN{{{args.n}}}"
    else:
        addition = f"N{{{args.n}}}X"

    count = modify_fasta(args.input, args.output, addition, args.position, args.replace_i)
    verb = "Prepended" if args.position == "start" else "Appended"
    print(f"Done. {verb} '{addition}' to {count} primer(s). Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
