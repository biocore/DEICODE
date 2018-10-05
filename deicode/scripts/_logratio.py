from biom import load_table
import pandas as pd
import numpy as np
import skbio
import os
from deicode.ratios import log_ratios
from deicode.preprocessing import rclr
import click
from gneiss.util import match
from skbio.stats.ordination import OrdinationResults

@click.command()
@click.option('--in_biom', help='Input table in biom format. (optional taxa in observation)')
@click.option('--in_map', help='Input metadata in qiime stype format.')
@click.option('--in_ord', help='RPCA output RPCA_Ordination.txt')
@click.option('--output_dir', help='Location of output files.')
@click.option('--axis', help='Axis of both ordinations to use default=0',default=0)
@click.option(
    '--min_sample_depth',
    default=500,
    help='Minimum Sample Sequencing Depth Cut Off default=500')
def logratio(in_biom: str, in_map: str,
             in_ord:str, output_dir: str,
             factor:str, min_sample_depth: int,
             axis:int) -> None:
    """ Runs log ratios on import features from RPCA output"""

    table = load_table(in_biom)
    metadata=pd.read_table(in_map)

    def sample_filter(val, id_, md): return sum(val) > min_sample_depth
    table = table.filter(sample_filter, axis='sample')
    table = table.to_dataframe().T.drop_duplicates()   
    # match tables
    table,metadata=match(table,metadata)
    # get loadings from ordination files
    sample_loading=OrdinationResults.read(in_ord).samples
    feature_loading=OrdinationResults.read(in_ord).biplot_scores

    try: 
        # try to inclide taxa if there 
        taxonomy=table.metadata_to_dataframe('observation') 
        logdf=log_ratios(table, metadata, feature_loading, 
                                 sample_loading, factor, taxa_tmp=taxonomy, 
                                 axis_sort=axis)
    except:
        # if not then just run with OTU ids
        logdf=log_ratios(table, metadata, feature_loading, 
                                 sample_loading, factor,
                                axis_sort=axis)

    logdf.to_csv(os.path.join(output_dir,'distance.txt'), sep='\t')

    return


if __name__ == '__main__':
    logratio()
