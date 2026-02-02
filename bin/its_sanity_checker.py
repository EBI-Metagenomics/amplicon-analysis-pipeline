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

TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD = 1000
MAPPING_PROPORTION_TEST_THRESHOLD = 0.50
RANK_PROPORTION_TEST_THRESHOLD = 0.50


def get_linecount(input_file: Path) -> int:
    """Get number of non-header lines in a file"""
    # the -1 is to not count headers
    linecount = sum(1 for _ in open(input_file, "r")) - 1
    return linecount


def tax_assignment_count_test(itsonedb_linecount: int, unite_linecount: int) -> bool:
    """
    Tax Assignment Count test will check whether the number of reads with assignments
    is above the `TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD` threshold
    """
    itsonedb_pass = itsonedb_linecount > TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD
    unite_pass = unite_linecount > TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD

    return itsonedb_pass or unite_pass


def mapping_proportion_test(
    itsonedb_linecount: int, unite_linecount: int, rrna_readcount: int
) -> tuple[bool, float, float]:
    """
    Proportion test will check whether the proportion of reads with assignments
    is above the `PROPORTION_TEST_THRESHOLD` threshold
    """

    # immediate failure of this test in these =0 scenarios
    if rrna_readcount == 0:
        return False, 0.0, 0.0
    if itsonedb_linecount == 0:
        logging.info("ITSoneDB Mapping Proportion: 0")
        itsonedb_pass = False
        itsonedb_reads_mapping = 0.0
    else:
        itsonedb_reads_mapping = itsonedb_linecount / float(rrna_readcount)
        logging.info(f"ITSoneDB Mapping Proportion: {itsonedb_reads_mapping}")
        itsonedb_pass = itsonedb_reads_mapping > MAPPING_PROPORTION_TEST_THRESHOLD
    if unite_linecount == 0:
        logging.info("UNITE Mapping Proportion: 0")
        unite_pass = False
        unite_reads_mapping = 0.0
    else:
        unite_reads_mapping = unite_linecount / float(rrna_readcount)
        logging.info(f"UNITE Mapping Proportion: {unite_reads_mapping}")
        unite_pass = unite_reads_mapping > MAPPING_PROPORTION_TEST_THRESHOLD

    return itsonedb_pass or unite_pass, itsonedb_reads_mapping, unite_reads_mapping


def rank_proportion_test(
    itsonedb_df: pd.DataFrame,
    unite_df: pd.DataFrame,
    itsonedb_linecount: int,
    unite_linecount: int,
) -> tuple[bool, float, float]:
    """
    Proportion test will check whether the proportion of reads with assignments
    below the rank of Kingdom is above the `RANK_PROPORTION_TEST_THRESHOLD` threshold
    """

    if itsonedb_linecount == 0:
        itsone_ranks_below_kingdom = 0
        logging.info("ITSoneDB Rank Proportion: 0")
    else:
        itsone_ranks_below_kingdom = len(
            itsonedb_df[itsonedb_df["taxon"].str.contains("p__")]
        ) / float(itsonedb_linecount)
        logging.info(f"ITSoneDB Rank Proportion: {itsone_ranks_below_kingdom}")

    if unite_linecount == 0:
        unite_ranks_below_kingdom = 0
        logging.info("UNITE Rank Proportion: 0")
    else:
        unite_ranks_below_kingdom = len(
            unite_df[unite_df["taxon"].str.contains("p__")]
        ) / float(unite_linecount)
        logging.info(f"UNITE Rank Proportion: {unite_ranks_below_kingdom}")

    itsonedb_pass = (
        itsone_ranks_below_kingdom > RANK_PROPORTION_TEST_THRESHOLD
        and itsonedb_linecount > TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD
    )
    unite_pass = (
        unite_ranks_below_kingdom > RANK_PROPORTION_TEST_THRESHOLD
        and unite_linecount > TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD
    )

    return (
        itsonedb_pass or unite_pass,
        itsone_ranks_below_kingdom,
        unite_ranks_below_kingdom,
    )


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
    Runs three sanity tests to verify whether reads are actually from ITS or from a different marker gene.
    Inputs are:
        - Uses mapping output files from MAPseq and UNITE + ITSoneDB reference databases
        - rRNA extraction output reads which are 'potential' ITS reads

    The tests are:
        - **Tax Assignment Count Test**: Check whether the number of reads with assignments
        is above the `TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD` threshold
        - **Mapping Proportion Test**: Check whether the proportion of reads with assignments
        is above the `PROPORTION_TEST_THRESHOLD` threshold
        - **Tax Assignment Count Test**: Check whether the proportion of reads with assignments
        below the rank of Kingdom is above the `RANK_PROPORTION_TEST_THRESHOLD` threshold

    Outputs are:
        - A `tsv` file containing the result for each test for the input
    """
    logging.info("Running ITS sanity checker on these inputs:")
    logging.info(
        f"{itsonedb_output=}, {unite_output=}, {rrna_extraction_output=}, {output_prefix=}"
    )

    results_dict = {}

    itsonedb_df = pd.read_csv(
        itsonedb_output, header=0, delim_whitespace=True, usecols=[12], names=["taxon"]
    )
    unite_df = pd.read_csv(
        unite_output, header=0, delim_whitespace=True, usecols=[12], names=["taxon"]
    )

    itsonedb_linecount = len(itsonedb_df)
    unite_linecount = len(unite_df)
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

    logging.info("Running Mapping Proportion Test")
    mapping_proportion_pass, itsonedb_reads_mapping, unite_reads_mapping = (
        mapping_proportion_test(itsonedb_linecount, unite_linecount, rrna_readcount)
    )
    results_dict["mapping_proportion_test"] = mapping_proportion_pass

    logging.info("Running Rank Proportion Test")

    rank_proportion_pass, itsone_ranks_below_kingdom, unite_ranks_below_kingdom = (
        rank_proportion_test(itsonedb_df, unite_df, itsonedb_linecount, unite_linecount)
    )
    results_dict["rank_proportion_test"] = rank_proportion_pass

    if mapping_proportion_pass:
        logging.info("Mapping Proportion test: PASSED")
    else:
        logging.info("Mapping Proportion test: FAILED")

    if tax_assignment_count_pass:
        logging.info("Tax Assignment Count test: PASSED")
    else:
        logging.info("Tax Assignment Count test: FAILED")

    if rank_proportion_pass:
        logging.info("Rank Proportion test: PASSED")
    else:
        logging.info("Rank Proportion test: FAILED")

    # Add a couple more numbers into the dictionary for debugging
    results_dict["rrna_readcount"] = rrna_readcount
    results_dict["itsonedb_linecount"] = itsonedb_linecount
    results_dict["unite_linecount"] = unite_linecount
    results_dict["itsonedb_reads_mapping"] = itsonedb_reads_mapping
    results_dict["unite_reads_mapping"] = unite_reads_mapping
    results_dict["itsone_ranks_below_kingdom"] = itsone_ranks_below_kingdom
    results_dict["unite_ranks_below_kingdom"] = unite_ranks_below_kingdom

    res_df = pd.DataFrame(results_dict, index=[0])

    out_path = f"{output_prefix}_its_sanity_check.json"
    logging.info(f"Saving results of tests to {out_path}")
    res_df.to_json(out_path, orient="records", double_precision=3, indent=2)


def main():
    its_sanity_checker()


if __name__ == "__main__":
    main()
