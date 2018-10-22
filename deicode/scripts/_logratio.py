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
@click.option('--in_ord', help='RPCA output RPCA_Ordination.txt')
@click.option('--output_dir', help='Location of output files.')
@click.option('--axis', help='Axis of both ordinations to use default=0',default=0)
@click.option('--n_lr', help='Number of log-ratios to comput default=10',default=10)
@click.option('--tax_level', help='If taxa included - choose level default=lowest',default=10)
def logratio(in_biom: str, in_ord:str, 
             output_dir: str, axis:int, 
             n_lr: int, tax_level: str) -> None:
    """ Runs log ratios on import features from RPCA output"""

    table = load_table(in_biom)
    tabledf = table.to_dataframe().T.drop_duplicates()
    # get loadings from ordination files
    sample_loading=OrdinationResults.read(in_ord).samples
    feature_loading=OrdinationResults.read(in_ord).features
    # match tables
    tabledf,feature_loading=match(tabledf.T,feature_loading)
    tabledf,sample_loading=match(tabledf.T,sample_loading)

    try:
        # try to inclide taxa if there
        taxonomy=table.metadata_to_dataframe('observation')
        logdf=log_ratios(tabledf, feature_loading, 
                                 sample_loading, taxa_tmp=taxonomy, 
                                 axis_sort=axis,N_show=n_lr, 
                                 level=tax_level)
    except:
        # if not then just run with OTU ids
        logdf=log_ratios(tabledf, feature_loading, 
                                sample_loading,
                                axis_sort=axis,N_show=n_lr,
                                level=tax_level)

    logdf.to_csv(os.path.join(output_dir,'Log_Ratios.txt'), sep='\t')

    return


if __name__ == '__main__':
    logratio()
