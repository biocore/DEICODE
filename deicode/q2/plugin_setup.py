# ----------------------------------------------------------------------------
# Copyright (c) 2016--, deicode development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin
import qiime2.sdk
from deicode import __version__
from deicode.rpca import rpca
from qiime2.plugin import (Properties, Int)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.ordination import PCoAResults

# citations = qiime2.plugin.Citations.load(
#             'citations.bib', package='deicode')

plugin = qiime2.plugin.Plugin(
    name='deicode',
    version=__version__,
    website="https://github.com/biocore/DEICODE",
    # citations=[citations['martino-unpublished']],
    short_description=('Plugin for Robust Aitchison PCA:'
                       'compositional biplots from sparse count data.'),
    description=('This is a QIIME 2 plugin supporting Robust Aitchison on '
                 'feature tables'),
    package='deicode')

plugin.methods.register_function(
    function=rpca,
    inputs={'table': FeatureTable[Frequency]},
    parameters={
        'rank': Int,
        'min_sample_count': Int,
        'min_feature_count': Int,
        'iterations': Int,
    },
    outputs=[
        ('biplot', PCoAResults % Properties("biplot")),
        ('distance_matrix', DistanceMatrix)
    ],
    input_descriptions={
        'table': 'Input table of counts.',
    },
    parameter_descriptions={
        'rank': ('The underlying low-rank'
                 ' structure (suggested: 1 < rank < 10)'),
        'min_sample_count': ('Minimum sum cutoff of'
                             ' sample across all features'),
        'min_feature_count': ('Minimum sum cutoff of'
                              ' features across all samples'),
        'iterations': ('The number of iterations to optimize the solution'
                       ' (suggested to be below 100; beware of overfitting)'),
    },
    output_descriptions={
        'biplot': ('A biplot of the (Robust Aitchison) RPCA feature loadings'),
        'distance_matrix': ('The Aitchison distance of'
                            'the sample loadings from RPCA.')
    },
    name='(Robust Aitchison) RPCA Biplot',
    description=("Performs robust center log-ratio transform "
                 "robust PCA and ranks the features by the "
                 "loadings of the resulting SVD."),
    citations=[]
)
