import os
import skbio
import click
import pandas as pd
from biom import load_table
from skbio import OrdinationResults
from deicode.rpca import rpca
from deicode.optspace import OptSpace
from deicode.preprocessing import rclr


@click.command()
@click.option('--in-biom', help='Input table in biom format.', required=True)
@click.option('--output-dir', help='Location of output files.', required=True)
@click.option(
    '--rank',
    default=3,
    help='Rank with which to run OptSpace [default: 3]')
@click.option(
    '--min-sample-count',
    default=500,
    help='Minimum sum cutoff of sample across all features [default: 500]')
@click.option(
    '--min-feature-count',
    default=10,
    help='Minimum sum cutoff of features across all samples [default: 10]')
@click.option(
    '--iterations',
    default=5,
    help='The number of iterations to optimize the solution'
         ' (suggested to be below 100; beware of overfitting)'
         ' [default: 5]')
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
    ord_res.write(os.path.join(output_dir, 'RPCA_Ordination.txt'))
    dist_res.write(os.path.join(output_dir, 'RPCA_distance.txt'))


if __name__ == '__main__':
    standalone_rpca()
