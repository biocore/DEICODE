# DEICODE changelog

## Version 0.2.4 (2020-01-07)

### Features

* Feature presence frequency filter [raised in issue #45](https://github.com/biocore/DEICODE/issues/45)

This new feature allows users to choose a number from 0-100, with decimals allowed. The filter
will remove features that are present in less than that percentage of samples. In practice
this filter is useful to filter out features that will be difficult to use for log-ratios
downstream in [Qurro](https://github.com/biocore/qurro).

* 'optspace' option in n-components to make an estimation on the rank parameter. See below in bug-fixes for more info.

### Bug fixes

* Axis order [raised in issue #52](https://github.com/biocore/DEICODE/issues/52) and [issue #32](https://github.com/biocore/DEICODE/issues/32)

This was partially fixed in PR #33 and PR #34 but there was a remaining bug
in the optspace function [here](https://github.com/biocore/DEICODE/blob/b1f37059b79f6ca2e2db11ba2fb7500c1a92f87e/deicode/optspace.py#L211-L212).
This bug was subtle and only caused the order to change periodically. Examples of this fix being tested on both
real and simulated data can be found [here](https://nbviewer.jupyter.org/github/cameronmartino/hartig-net/blob/master/percent-explained/real-data-fraction-var.ipynb) and [here](https://nbviewer.jupyter.org/github/cameronmartino/hartig-net/blob/master/percent-explained/simulation-fraction-var.ipynb) respectively. 

* Fraction total variance explained [raised in issue #53](https://github.com/biocore/DEICODE/issues/53)

The fraction variance explained currently is calculated by the sum of squares based on the number of singular values. The number of singular values changes based on the rank input. This causes the resulting fraction variance explained to change dramatically as the rank is increased or decreased. This has been brought up multiple times on the QIIME2 forum. For example, see [here](). Additionally, the correct rank that actually explains the total variance is not clear (i.e. What rank should I choose?). To solve both of these problems a rank estimation is made based on part C of Keshavan, R. H., Montanari, A. & Oh, S. (Low-rank matrix completion with noisy observations: A quantitative comparison. in 2009 47th Annual Allerton Conference on Communication, Control, and Computing (Allerton) 1216–1222 (2009)) is used. This can can be enabled by setting n-components to 'optspace'. To prevent user confusion the default was not changed in this version. However, a future warning was added to warn users that in the next version 'optspace' based rank estimation will be the default.

## Version 0.2.3 (2019-6-18)

### Backward-incompatible changes [stable]

* The following option names have changed when running DEICODE:
  (this is true for both QIIME2, python, and the standalone interfaces)

| Original Name        | New Name             |
| -------------------- | -------------------- |
| `--rank`             | `--n_components`     |
| `--iterations`       | `--max_iterations`   |

* `deicode.optspace` has been replaced by `deicode.matrix_completion`.
  * `deicode.optspace.OptSapce` has been replaced by `deicode.matrix_completion.MatrixCompletion`.
* `deicode._optspace` has been replaced by `deicode.optspace`.
  * `deicode._optspace.optspace` has been replaced by `deicode.optspace.OptSpace`.

### Performance enhancements

* Updated `deicode.optspace.OptSpace.svd_sort` to be readable and included
  comments to make the code easier to follow
* Added several tests that were missing for ValueError(s) to get better code coverage

### Bug fixes

* Arrow scaling issues [first raised here](https://forum.qiime2.org/t/deicode-rpca-biplot-are-the-vectors-distorting-the-plot/9497).
  * this was caused by some bit of code that was borrowed from sklearn in
    `deicode.optspace.OptSpace.svd_sort`. Just deleting a few lines fixed it.

## Version 0.2.2 (2019-4-22)

### Features

* Ensure sorting that the eigenvalues is being
  is from largest to smallest. This is done for
  the singular value matrix, given as the variable
  s. This involves the sorting s by the diagonal
  while also ordering off diagonal elements. Additionally,
  a new function was added to ensure that the sorted
  eigenvalues are also ordered in the U and V
  loading values in the SVD. This function also
  ensures that after sorting, the values of the SVD
  are deterministic. This is all implemented in the
  function located in _optspace.py called svd_sort.
  This methodology is also performed in the code
  for PCA in scikit learn. To do this I pulled
  from this code and noted this in the code
  comments with direct line links in scikit-learn. 

* Tests for the function svd_sort described above
  were added. The test was added under the location
  tests/test_optspace.py in the function given by
  the function test_optspace_ordering.

## Version 0.2.1 (2019-4-15)

Implemented in [PR#33](https://github.com/biocore/DEICODE/pull/33).

### Bug fixes

* Eigenvalues were out of order, single line problem in optspace.py

## Version 0.2.0 (2019-4-8)

Implemented in [PR#29](https://github.com/biocore/DEICODE/pull/29).

### Features

* The `--min-feature-count` and `--iterations` options are now usable when
  running DEICODE outside of QIIME 2.

* When running DEICODE outside of QIIME 2, it will no longer just raise an error if
  the specified output directory doesn't exist (it'll automatically try to
  create the output directory, even creating multiple levels of directories if
  needed).

* Options' default values (if applicable) are now shown when running
  `deicode --help`.

* The `--in-biom` and `--output-dir` options are now marked as required when
  running DEICODE outside of QIIME 2 (so just running, e.g., `deicode` without
  any options will give you a clearer error message than before).

* The RPCA functionality is now exposed via the `deicode.rpca` module,
  which contains an `rpca()` function.

* Duplicate indices and columns will cause a ValueError. Previously
  the script `deicode.scripts.rpca.py` would just drop any duplicates.

### Backward-incompatible changes [stable]

* The following option names have changed when running DEICODE outside of QIIME
  2:

| Original Name        | New Name             |
| -------------------- | -------------------- |
| `--in_biom`          | `--in-biom`          |
| `--output_dir`       | `--output-dir`       |
| `--min_sample_depth` | `--min-sample-count` |

* `deicode.scripts._rpca` has been replaced by `deicode.scripts._standalone_rpca`.
  * Similarly, the `rpca()` function within `deicode.scripts._rpca` has been
    replaced by the `standalone_rpca()` function in
    `deicode.scripts._standalone_rpca`.
  * `deicode.preprocessing.inverse_rclr` was removed along with its
     tests. This code was redundant with `skbio.stats.composition.clr_inv`,
     which can be performed on clr-transformed data. Furthermore,
     this inverse is a holdover from old versions of DEICODE
     where we directly interpreted the imputation and is no
     longer useful for the output. 

* `deicode.ratios.py` was removed. This (untested) code was an unfinished
   feature that will be replaced by rankratioviz. This code was only used
   in the visualizations used in the manuscript and will be stored in the
   DEICODE-benchmarking repository where it is actually used.

### Backward-incompatible changes [experimental]

### Performance enhancements

* Removed scikit-learn dependency.

### Bug fixes

* Some of the QIIME 2 RPCA behavior was not mirrored perfectly in the non-QIIME
  2 RPCA code. Here is a list of the "new" things done by the non-Q2 RPCA code
  that it didn't do before:
  * Uses `--min-feature-count` with a default value of `10`. Previously, the
    non-Q2 RPCA code didn't do this filtering step at all.
  * Adds a PC3 containing zeros if the `rank` is set to `2` (to support
    visualizing these biplots in Emperor).

* A minimum value of `2` is now enforced for the `--rank` option.

* A minimum value of `1` is now enforced for the `--iterations` option.

* Fixed the test in `deicode/scripts/tests/` to check the correct output
  files produced by DEICODE (previously, this test was looking at the incorrect
  files).

* Fixed a test in `deicode/q2/tests/` to correctly check for NaNs in the
  ordination produced by DEICODE (previously, this test was using python's
  built-in `any()` function instead of pandas' `.any()` function, which
  resulted in the test being incorrect).

* Iteration in `deicode/_optspace.py` indexing was off see @fedarko's 
  comment in PR #29. This causes the iteration to be one less than the
  input, this should not have had an impact any results.

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

* Since `deicode.rpca` is now used by both the QIIME 2 and non-QIIME 2 code,
  the amount of redundant code has decreased. This should simplify DEICODE's
  codebase.

* Shared RPCA parameter settings between the QIIME 2 and non-QIIME 2 code
  (descriptions and default values) are now stored in `deicode._rpca_defaults`.
  This further cuts down on redundancy: developers now only have to update the
  description or default value of an option in this one place.

* The files produced by the non-QIIME 2 code have been renamed as follows to be
  consistent with the data inside the artifacts produced by the QIIME 2 code:

| Original Name         | New Name              |
| --------------------- | --------------------- |
| `RPCA_Ordination.txt` | `ordination.txt`      |
| `RPCA_distance.txt`   | `distance-matrix.tsv` |

* The `skbio.OrdinationResults` short and long method names for ordinations
  produced by the non-QIIME 2 code have changed as follows, in order to be
  consistent with the ordinations produced by the QIIME 2 code:

| Original Name                   | New Name                         |
| ------------------------------- | -------------------------------- |
| `PCoA`                          | `rpca_biplot`                    |
| `Principal Coordinate Analysis` | `(Robust Aitchison) RPCA Biplot` |

  * These changes shouldn't actually impact anything within the
    ordination files, since as of writing the [scikit-bio ordination format](http://scikit-bio.org/docs/latest/generated/skbio.io.format.ordination.html) doesn't include either of the method names. Just listing this here to be safe.

* Various typo fixes

* Various RPCA test code enhancements

* A citation of DEICODE's
  [published paper](https://msystems.asm.org/content/4/1/e00016-19) is now
  included in the citations of QIIME 2 artifacts generated with DEICODE.

## Version 0.1.6 (2019-3-8)

In PR#27

### Features

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

### Bug fixes

* Centered the feature and sample loadings for biplot visualization (issue #26).

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

* Fixed the broken image links.

## Version 0.1.5

Original "working" code

### Features

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

### Bug fixes

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous
