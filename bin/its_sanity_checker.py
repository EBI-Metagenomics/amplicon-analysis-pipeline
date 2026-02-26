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

import json
import logging

import click
import pandas as pd
from Bio import SeqIO

logging.basicConfig(level=logging.DEBUG)

TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD = 1000
MAPPING_PROPORTION_TEST_THRESHOLD = 0.50
RANK_PROPORTION_TEST_THRESHOLD = 0.50


def count_sequences_in_file(filepath: str) -> int:
    """
    Count the number of sequences in a FASTA or FASTQ file using SeqIO.
    """
    try:
        # Try FASTA first
        count = sum(1 for _ in SeqIO.parse(filepath, "fasta"))
        if count > 0:
            return count
    except Exception:
        pass

    try:
        # Try FASTQ
        count = sum(1 for _ in SeqIO.parse(filepath, "fastq"))
        if count > 0:
            return count
    except Exception:
        pass

    # If both fail, raise an error
    raise ValueError(f"Could not parse {filepath} as FASTA or FASTQ format")


def read_mseq_file(filepath: str) -> pd.DataFrame:
    """
    Read a MAPseq output file into a DataFrame with a 'taxon' column.
    Returns an empty DataFrame if the file is empty or cannot be parsed.
    """
    try:
        return pd.read_csv(
            filepath,
            sep=r"\s+",
            comment="#",
            header=None,
            usecols=[12],
            names=["taxon"],
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["taxon"])


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


def _run_sanity_tests(filepaths: dict[str, str]) -> dict:
    """Run all three sanity tests and return the results dictionary."""
    # Read MAPseq output files
    logging.info("Reading MAPseq assignment files...")
    itsonedb_df = read_mseq_file(filepaths["ITSoneDB"])
    itsonedb_linecount = len(itsonedb_df)
    logging.info("ITSoneDB assignment count: %d", itsonedb_linecount)

    unite_df = read_mseq_file(filepaths["UNITE"])
    unite_linecount = len(unite_df)
    logging.info("UNITE assignment count: %d", unite_linecount)

    try:
        rrna_readcount = count_sequences_in_file(filepaths["Rfam_SSU_LSU"])
        logging.info("Potential ITS reads count: %d", rrna_readcount)
    except Exception as e:
        logging.error("Failed to count Rfam SSU+LSU sequences: %s", e)
        rrna_readcount = 0

    results_dict = {}

    # Tax Assignment Count Test
    logging.info("Running Tax Assignment Count Test")
    tax_assignment_count_pass = tax_assignment_count_test(
        itsonedb_linecount, unite_linecount
    )
    results_dict["tax_assignment_count_test"] = tax_assignment_count_pass

    # Mapping Proportion Test
    logging.info("Running Mapping Proportion Test")
    mapping_proportion_pass, itsonedb_reads_mapping, unite_reads_mapping = (
        mapping_proportion_test(itsonedb_linecount, unite_linecount, rrna_readcount)
    )
    results_dict["mapping_proportion_test"] = mapping_proportion_pass

    # Rank Proportion Test
    logging.info("Running Rank Proportion Test")
    rank_proportion_pass, itsone_ranks_below_kingdom, unite_ranks_below_kingdom = (
        rank_proportion_test(
            itsonedb_df, unite_df, itsonedb_linecount, unite_linecount
        )
    )
    results_dict["rank_proportion_test"] = rank_proportion_pass

    # Log test results
    for name, passed in [
        ("Mapping Proportion", mapping_proportion_pass),
        ("Tax Assignment Count", tax_assignment_count_pass),
        ("Rank Proportion", rank_proportion_pass),
    ]:
        logging.info("%s test: %s", name, "PASSED" if passed else "FAILED")

    # Add counts into the dictionary for debugging
    results_dict["rrna_readcount"] = rrna_readcount
    results_dict["itsonedb_linecount"] = itsonedb_linecount
    results_dict["unite_linecount"] = unite_linecount
    results_dict["itsonedb_reads_mapping"] = itsonedb_reads_mapping
    results_dict["unite_reads_mapping"] = unite_reads_mapping
    results_dict["itsone_ranks_below_kingdom"] = itsone_ranks_below_kingdom
    results_dict["unite_ranks_below_kingdom"] = unite_ranks_below_kingdom

    return results_dict


############## ITS Sanity Checker ##############
@click.command(
    short_help="Sanity check whether a potential ITS result isn't just a different marker gene",
)
@click.option(
    "--read_assignments",
    required=True,
    help="JSON file containing file paths to read assignment results keyed by database name",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "-p",
    "--output_prefix",
    required=True,
    help="Prefix to sanity checker output",
    type=str,
)
def its_sanity_checker(
    read_assignments: str,
    output_prefix: str,
) -> None:
    """
    Runs three sanity tests to verify whether reads are actually from ITS or from a different marker gene.
    Inputs are:
        - A JSON file containing file paths to:
          - ITSone: MAPseq output file for ITSoneDB
          - UNITE: MAPseq output file for UNITE
          - Rfam_SSU_LSU: FASTA or FASTQ file containing Rfam-classified SSU+LSU reads

    The tests are:
        - **Tax Assignment Count Test**: Check whether the number of reads with assignments
        is above the `TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD` threshold
        - **Mapping Proportion Test**: Check whether the proportion of reads with assignments
        is above the `PROPORTION_TEST_THRESHOLD` threshold
        - **Rank Proportion Test**: Check whether the proportion of reads with assignments
        below the rank of Kingdom is above the `RANK_PROPORTION_TEST_THRESHOLD` threshold

    Outputs are:
        - A `json` file containing the result for each test for the input
        - A `tsv` file containing the result for each test for the input. Mainly used for multiqc in the amplicon analysis pipeline
    """
    logging.info("Running ITS sanity checker on these inputs:")
    logging.info(f"{read_assignments=}, {output_prefix=}")

    with open(read_assignments) as f:
        filepaths = json.load(f)

    # Check for required database keys
    db_keys = ["ITSoneDB", "UNITE"]
    missing_db_keys = [k for k in db_keys if k not in filepaths]
    if missing_db_keys:
        logging.warning(
            "Missing database results in read assignments JSON: %s. "
            "Outputting negative sanity check result.",
            missing_db_keys,
        )
        results_dict = {
            "tax_assignment_count_test": False,
            "mapping_proportion_test": False,
            "rank_proportion_test": False,
            "rrna_readcount": 0,
            "itsonedb_linecount": 0,
            "unite_linecount": 0,
            "itsonedb_reads_mapping": 0.0,
            "unite_reads_mapping": 0.0,
            "itsone_ranks_below_kingdom": 0.0,
            "unite_ranks_below_kingdom": 0.0,
        }
    else:
        results_dict = _run_sanity_tests(filepaths)

    # set the prefix as the index, implied that it's the run ID for our purposes
    results_dict["run"] = output_prefix

    out_path = f"{output_prefix}_its_sanity_check.json"
    # multiqc processing is easier if the file is a tsv/csv
    out_path_mqc = f"{output_prefix}_its_sanity_check_mqc.tsv"

    logging.info(f"Saving results of tests to {out_path} and {out_path_mqc}")

    # Write JSON output
    with open(out_path, "w") as f:
        json.dump([results_dict], f, indent=2)

    # Write TSV output for multiqc
    header = "\t".join(k for k in results_dict if k != "run")
    values = "\t".join(str(results_dict[k]) for k in results_dict if k != "run")
    with open(out_path_mqc, "w") as f:
        f.write(f"run\t{header}\n")
        f.write(f"{output_prefix}\t{values}\n")


def main():
    its_sanity_checker()


if __name__ == "__main__":
    main()
