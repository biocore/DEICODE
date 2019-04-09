# Configuration file where you can set the parameter default values and
# descriptions. This is used by both the standalone RPCA and QIIME 2 RPCA sides
# of DEICODE.
DEFAULT_RANK = 3
DEFAULT_MSC = 500
DEFAULT_MFC = 10
DEFAULT_ITERATIONS = 5

DESC_RANK = ("The underlying low-rank structure (suggested: 1 < rank < 10)"
             " [minimum 2]")
DESC_MSC = "Minimum sum cutoff of sample across all features"
DESC_MFC = "Minimum sum cutoff of features across all samples"
DESC_ITERATIONS = ("The number of iterations to optimize the solution"
                   " (suggested to be below 100; beware of overfitting)"
                   " [minimum 1]")
