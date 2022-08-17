import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


def biplot(axis1, axis2, sload, fload, hue, ax,
           n_arrow= 15, level=2, cmap = 'Set1'):
    
    """
    This function is a helper for the tutorial. It
    takes ordination inputs from the tutorial and
    plots a compositional biplot.
    
    Parameters
    ----------
    axis1: str - first x-PC axis

    axis2: str - second y-PC axis
    
    sload: pd.DataFrame - sample ordination

    fload: pd.DataFrame - feature ordination

    hue: str - subject groupings
        
    ax: matplitlib subplot

    n_arrow: int - number of arrows to plot

    level: int - taxonomic level to collapse

    Returns
    ----------
    ax: matplitlib subplot

    """

    # sort the arrows by thier PC1, PC2 mag.
    fload['magnitude'] = np.linalg.norm(fload[[axis1, axis2]], axis=1)
    fload = fload.sort_values('magnitude', ascending=False)
    
    # get a taxonomic level split
    def _collapse(tax, max_observed_level = 7):
        tax = [x.strip() for x in tax.split(';')]
        if len(tax) < max_observed_level:
            padding = ['__'] * (max_observed_level - len(tax))
            tax.extend(padding)
        return ';'.join(tax[:level])
    fload['level'] = fload.Taxon.apply(_collapse)
    
    # get cmap for arrows
    arrow_cmap = set(fload.loc[fload.index[: n_arrow], 'level'].values)
    colors_ = plt.get_cmap(cmap, len(arrow_cmap))
    arrow_cmap = {v:colors_(i) for i, v in enumerate(arrow_cmap)}

    # plot the sample loadings
    sns.scatterplot(axis1, axis2,
                    data = sload,
                    hue = hue,
                    palette=cmap,
                    ax=ax)
    # plot the arrows
    legend_arrows = {}
    for i in range(n_arrow):
        # add arrow
        ind_ = fload.index[i]
        arrow_  = ax.arrow(0, 0, fload.loc[ind_, axis1],
                        fload.loc[ind_, axis2],
                        color=arrow_cmap[fload.loc[ind_, 'level']],
                        width=0.001, head_width=0.005)
        legend_arrows[fload.loc[ind_, 'level']] = arrow_
    
    # add legend
    leg1 = ax.legend(loc='center left',
                     bbox_to_anchor=(1, .8),
                     title="sample legend (dots)")
    ax.legend(list(legend_arrows.values()),
              list(legend_arrows.keys()), 
              loc='center left',
              bbox_to_anchor=(1, 0.5),
              title="feature legend (arrows)")
    ax.add_artist(leg1)
        
    return ax
