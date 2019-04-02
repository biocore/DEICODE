import os
import click
from biom import load_table
from deicode.rpca import rpca
from deicode._rpca_defaults import (DEFAULT_RANK, DEFAULT_MSC, DEFAULT_MFC,
                                    DEFAULT_ITERATIONS, DESC_RANK, DESC_MSC,
                                    DESC_MFC, DESC_ITERATIONS)


@click.command()
@click.option('--in-biom', help='Input table in biom format.', required=True)
@click.option('--output-dir', help='Location of output files.', required=True)
@click.option(
    '--rank',
    default=DEFAULT_RANK,
    show_default=True,
    help=DESC_RANK)
@click.option(
    '--min-sample-count',
    default=DEFAULT_MSC,
    show_default=True,
    help=DESC_MSC)
@click.option(
    '--min-feature-count',
    default=DEFAULT_MFC,
    show_default=True,
    help=DESC_MFC)
@click.option(
    '--iterations',
    default=DEFAULT_ITERATIONS,
    show_default=True,
    help=DESC_ITERATIONS)
def standalone_rpca(in_biom: str, output_dir: str, rank: int,
                    min_sample_count: int, min_feature_count: int,
                    iterations: int) -> None:
    """Runs RPCA with an rclr preprocessing step."""

    # import table
    table = load_table(in_biom)

    ord_res, dist_res = rpca(table, rank, min_sample_count, min_feature_count,
                             iterations)

    # If it doesn't already exist, create the output directory.
    # Note that there is technically a race condition here: it's ostensibly
    # possible that some process could delete the output directory after we
    # check that it exists here but before we write the output files to it.
    # However, in this case, we'd just get an error from skbio.io.util.open()
    # (which is called by skbio.OrdinationResults.write()), which makes sense.
    os.makedirs(output_dir, exist_ok=True)

    # write files to output directory
    # Note that this will overwrite files in the output directory that share
    # these filenames (analogous to QIIME 2's behavior if you specify the
    # --o-biplot and --o-distance-matrix options, but differing from QIIME 2's
    # behavior if you specify --output-dir instead).
    ord_res.write(os.path.join(output_dir, 'ordination.txt'))
    dist_res.write(os.path.join(output_dir, 'distance-matrix.tsv'))


if __name__ == '__main__':
    standalone_rpca()
