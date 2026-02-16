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

logging.basicConfig(level=logging.DEBUG)

TAX_ASSIGNMENT_COUNT_TEST_THRESHOLD = 1000
MAPPING_PROPORTION_TEST_THRESHOLD = 0.50
RANK_PROPORTION_TEST_THRESHOLD = 0.50



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
    short_help="Sanity check whether a potential ITS result isn't just a different marker gene",
)
@click.option(
    "--read_assignments",
    required=True,
    help="JSON file containing read assignment counts and file paths keyed by database name",
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
        - A JSON file containing read assignment counts (ITSone, UNITE, Rfam_SSU_LSU)
          and MAPseq output file paths (ITSone_fp, UNITE_fp) keyed by database name

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
    logging.info(
        f"{read_assignments=}, {output_prefix=}"
    )

    with open(read_assignments) as f:
        counts = json.load(f)

    # Log any missing keys
    expected_keys = ["ITSone", "UNITE", "Rfam_SSU_LSU", "ITSone_fp", "UNITE_fp"]
    missing_keys = [k for k in expected_keys if k not in counts]
    if missing_keys:
        logging.warning(f"Missing keys in read assignments JSON: {missing_keys}")

    itsonedb_linecount = counts.get("ITSone", 0)
    unite_linecount = counts.get("UNITE", 0)
    rrna_readcount = counts.get("Rfam_SSU_LSU", 0)

    logging.info(f"Potential ITS reads count: {rrna_readcount}")
    logging.info(f"ITSoneDB assignment count: {itsonedb_linecount}")
    logging.info(f"UNITE assignment count: {unite_linecount}")

    results_dict = {}

    # Tax Assignment Count Test requires: ITSone, UNITE
    logging.info("Running Tax Assignment Count Test")
    if "ITSone" in counts and "UNITE" in counts:
        tax_assignment_count_pass = tax_assignment_count_test(
            itsonedb_linecount, unite_linecount
        )
    else:
        logging.warning("Tax Assignment Count Test FAILED: missing required keys "
                        f"{[k for k in ['ITSone', 'UNITE'] if k not in counts]}")
        tax_assignment_count_pass = False
    results_dict["tax_assignment_count_test"] = tax_assignment_count_pass

    # Mapping Proportion Test requires: ITSone, UNITE, Rfam_SSU_LSU
    logging.info("Running Mapping Proportion Test")
    mapping_proportion_keys = ["ITSone", "UNITE", "Rfam_SSU_LSU"]
    if all(k in counts for k in mapping_proportion_keys):
        mapping_proportion_pass, itsonedb_reads_mapping, unite_reads_mapping = (
            mapping_proportion_test(itsonedb_linecount, unite_linecount, rrna_readcount)
        )
    else:
        logging.warning("Mapping Proportion Test FAILED: missing required keys "
                        f"{[k for k in mapping_proportion_keys if k not in counts]}")
        mapping_proportion_pass = False
        itsonedb_reads_mapping = 0.0
        unite_reads_mapping = 0.0
    results_dict["mapping_proportion_test"] = mapping_proportion_pass

    # Rank Proportion Test requires: ITSone, UNITE, ITSone_fp, UNITE_fp
    logging.info("Running Rank Proportion Test")
    rank_proportion_keys = ["ITSone", "UNITE", "ITSone_fp", "UNITE_fp"]
    if all(k in counts for k in rank_proportion_keys):
        itsonedb_df = pd.read_csv(
            counts["ITSone_fp"], header=0, delim_whitespace=True, usecols=[12], names=["taxon"]
        )
        unite_df = pd.read_csv(
            counts["UNITE_fp"], header=0, delim_whitespace=True, usecols=[12], names=["taxon"]
        )
        rank_proportion_pass, itsone_ranks_below_kingdom, unite_ranks_below_kingdom = (
            rank_proportion_test(itsonedb_df, unite_df, itsonedb_linecount, unite_linecount)
        )
    else:
        logging.warning("Rank Proportion Test FAILED: missing required keys "
                        f"{[k for k in rank_proportion_keys if k not in counts]}")
        rank_proportion_pass = False
        itsone_ranks_below_kingdom = 0.0
        unite_ranks_below_kingdom = 0.0
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
