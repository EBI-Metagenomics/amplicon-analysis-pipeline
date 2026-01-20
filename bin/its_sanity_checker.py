#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
from pathlib import Path

from Bio import SeqIO
import click
import pandas as pd

logging.basicConfig(level=logging.DEBUG)

TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD = 500
# TODO: increase this threshold, small for testing atm
PROPORTION_TEST_THRESHOLD = 0.10


def get_linecount(input_file: Path) -> int:
    """Get number of lines in a file"""
    linecount = sum(1 for _ in open(input_file, "r"))
    return linecount


def tax_assignment_count_test(itsonedb_linecount: int, unite_linecount: int) -> bool:
    """
    Tax Assignment Count test will check whether the number of reads with assignments
    is above the `TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD` threshold
    """
    itsonedb_pass = itsonedb_linecount > TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD
    unite_pass = unite_linecount > TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD

    # both DBs have to be above threshold to pass
    if itsonedb_pass or unite_pass:
        return True
    else:
        return False


def proportion_test(
    itsonedb_linecount: int, unite_linecount: int, rrna_readcount: int
) -> bool:
    """
    Proportion test will check whether the proportion of reads with assignments
    is above the `PROPORTION_TEST_THRESHOLD` threshold
    """

    # immediate failure of this test in these =0 scenarios
    if rrna_readcount == 0:
        return False
    if itsonedb_linecount == 0 and unite_linecount == 0:
        return False

    itsonedb_pass = (
        itsonedb_linecount / float(rrna_readcount) > PROPORTION_TEST_THRESHOLD
    )
    unite_pass = unite_linecount / float(rrna_readcount) > PROPORTION_TEST_THRESHOLD

    # both DBs have to be above threshold to pass
    if itsonedb_pass or unite_pass:
        return True
    else:
        return False


############## ITS Sanity Checker ##############
@click.command(
    options_metavar="--itsonedb_output <itsonedb_output> --unite_output <itsonedb_output> -r <reads> -p <output_prefix>",
    short_help="Sanity check whether a potential ITS result isn't just a different marker gene",
)
@click.option(
    "--itsonedb_output",
    required=True,
    help="Output from MAPseq containing ITSoneDB taxonomic assignments",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "--unite_output",
    required=True,
    help="Output from MAPseq containing UNITE taxonomic assignments",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-r",
    "--rrna_extraction_output",
    required=True,
    help="Output reads FASTA file from `MASK_FASTA_SWF:FILTER_MASKED_N` module which contains potential ITS assignments",
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
)
@click.option(
    "-p",
    "--output_prefix",
    required=True,
    help="Prefix to sanity checker output",
    type=str,
)
def its_sanity_checker(
    itsonedb_output: Path,
    unite_output: Path,
    rrna_extraction_output: Path,
    output_prefix: str,
) -> None:
    """
    Runs some sanity tests to verify whether a particular run is actually ITS
    or just a different marker gene
    """
    logging.info("Running ITS sanity checker on these inputs:")
    logging.info(
        f"{itsonedb_output=}, {unite_output=}, {rrna_extraction_output=}, {output_prefix=}"
    )

    results_dict = {}

    itsonedb_linecount = get_linecount(itsonedb_output)
    unite_linecount = get_linecount(unite_output)
    rrna_readcount = 0

    # get number of reads rather than number of lines
    with open(rrna_extraction_output) as fr:
        for _ in SeqIO.parse(fr, "fasta"):
            rrna_readcount += 1

    logging.info(f"Potential ITS reads count: {rrna_readcount}")
    logging.info(f"ITSoneDB assignment count: {itsonedb_linecount}")
    logging.info(f"UNITE assignment count: {unite_linecount}")

    logging.info("Running Tax Assignment Count Test")
    tax_assignment_count_pass = tax_assignment_count_test(
        itsonedb_linecount, unite_linecount
    )
    results_dict["tax_assignment_count_test"] = tax_assignment_count_pass

    logging.info("Running Proportion Test")
    proportion_pass = proportion_test(
        itsonedb_linecount, unite_linecount, rrna_readcount
    )

    results_dict["proportion_test"] = proportion_pass
    if proportion_pass:
        logging.info("Proportion test: PASSED")
    else:
        logging.info("Proportion test: FAILED")

    if tax_assignment_count_pass:
        logging.info("Tax Assignment Count test: PASSED")
    else:
        logging.info("Tax Assignment Count test: FAILED")

    res_df = pd.DataFrame(results_dict, index=[0])

    out_path = f"{output_prefix}_its_sanity_check.tsv"
    logging.info(f"Saving results of tests to {out_path}")
    res_df.to_csv(out_path, sep="\t", index=False)


def main():
    its_sanity_checker()


if __name__ == "__main__":
    main()
