# DEICODE changelog

## Version 0.1.7 (2019-4-2)

Implemented in [PR#29](https://github.com/biocore/DEICODE/pull/29).

### Features

* The `--min-feature-count` and `--iterations` options are now usable when
  running DEICODE outside of QIIME 2.

* When running DEICODE outside of QIIME 2, it will no longer raise an error if
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

### Backward-incompatible changes [stable]

* The following option names have changed when running DEICODE outside of QIIME
  2:
  | Original Name        | New Name             |
  |----------------------|----------------------|
  | `--in_biom`          | `--in-biom`          |
  | `--output_dir`       | `--output-dir`       |
  | `--min_sample_depth` | `--min-sample-count` |

* `deicode.scripts._rpca` has been replaced by `deicode.scripts._standalone_rpca`.
  * Similarly, the `rpca()` function within `deicode.scripts._rpca` has been
    replaced by the `standalone_rpca()` function in
    `deicode.scripts._standalone_rpca`.

### Backward-incompatible changes [experimental]

### Performance enhancements

### Bug fixes

* Some of the QIIME 2 RPCA behavior was not mirrored perfectly in the non-QIIME
  2 RPCA code. So there are a few "new" things done by the non-Q2 RPCA code that
  it didn't do before:
  * Uses `--min-feature-count` with a default value of `10`. Previously, the
    non-Q2 RPCA code didn't do this filtering step at all.
  * Adds a PC3 containing zeros if the `rank` is set to `2`.

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

* Since `deicode.rpca` is now used by both the QIIME 2 and non-QIIME 2 code,
  the amount of redundant code has decreased. The primary benefit of this is
  that it should maintaining both ways of running DEICODE easier.

* Shared RPCA parameter settings between the QIIME 2 and non-QIIME 2 code
  (descriptions and default values) are now stored in `deicode._rpca_defaults`.
  This further cuts down on redundancy: developers now only have to update the
  description or default value of an option in this one place.

* The `skbio.OrdinationResults` short and long method names for ordinations
  produced by the non-QIIME 2 code have changed as follows, in order to be
  consistent with the ordinations produced by the QIIME 2 code:
  | Original Name                   | New Name                         |
  |---------------------------------|----------------------------------|
  | `PCoA`                          | `rpca_biplot`                    |
  | `Principal Coordinate Analysis` | `(Robust Aitchison) RPCA Biplot` |
  * These changes shouldn't actually impact anything within the
    `RPCA_Ordination.txt` files, since as of writing the [scikit-bio ordination format](http://scikit-bio.org/docs/latest/generated/skbio.io.format.ordination.html) doesn't include either of the method names. Just listing this here to be safe.

* Various typo fixes

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