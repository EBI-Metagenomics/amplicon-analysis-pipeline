import logging
from pathlib import Path


import click
from Bio import SeqIO

logging.basicConfig(level=logging.DEBUG)


@click.command(
    options_metavar="-m <mapseq_output> -r <reads> -p <output_prefix>",
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
    Runs a bunch of sanity tests to verify whether a particular run is actually ITS
    or just a different marker gene
    # TODO: do these docstring params
    :param mapseq_output: Description
    :type mapseq_output: Path
    :param rrna_extraction_output: Description
    :type rrna_extraction_output: Path
    :param output_prefix: Description
    :type output_prefix: str
    """
    logging.info("Running ITS sanity checker on these inputs:")
    logging.info(
        f"{itsonedb_output=}, {unite_output=}, {rrna_extraction_output=}, {output_prefix=}"
    )

    itsonedb_linecount = get_linecount(itsonedb_output)
    unite_linecount = get_linecount(unite_output)
    rrna_linecount = 0

    with open(rrna_extraction_output) as fr:
        for _ in SeqIO.parse(fr, "fasta"):
            rrna_linecount += 1

    logging.info(f"{rrna_linecount}, {itsonedb_linecount}, {unite_linecount}")


def get_linecount(input_file: Path) -> int:
    linecount = sum(1 for i in open(input_file, "r"))
    return linecount


def main():
    its_sanity_checker()


if __name__ == "__main__":
    main()
