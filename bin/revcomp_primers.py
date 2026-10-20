#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

from Bio.Seq import Seq
from Bio import SeqIO
import argparse


def reverse_complement_iupac(sequence):
    return str(Seq(sequence.strip().upper()).reverse_complement())


def main():
    parser = argparse.ArgumentParser(
        description="Reverse complement sequences in a FASTA file while preserving IDs"
    )

    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output FASTA file")

    args = parser.parse_args()

    records = SeqIO.parse(args.input, "fasta")

    with open(args.output, "w") as output_handle:
        for record in records:
            # Keep the FASTA ID/header unchanged
            # and reverse-complement only the sequence
            record.seq = record.seq.reverse_complement()

            SeqIO.write(record, output_handle, "fasta")

    print(f"Done! Written to {args.output}")


if __name__ == "__main__":
    main()

