from __future__ import division

try:
    import numpy as np
except ImportError:
    print('Unable to import numpy.')
try:
    from sklearn.preprocessing import Imputer
except ImportError:
    print('Unable to import Imputer from sklearn.preprocessing.')

from sklearn.preprocessing import Imputer
import numpy as np

class base(object):
    
    '''''
        This contains base bench marks for impuations (all methods should far exceed the accuracy of theese)
        '''''
    
    #Zero fill
    def zeros(A):
        return np.nan_to_num(A)
    #return A
    
    #mean fill benchmark
    def mean_fill(A):
        imp =Imputer(missing_values='NaN', strategy='mean', axis=0, verbose=0, copy=True)
        imp.fit(A)
        return imp.transform(A)
    
    #median fill benchmark
    def median_fill(A):
        imp =Imputer(missing_values='NaN', strategy='median', axis=0, verbose=0, copy=True)
        imp.fit(A)
        return imp.transform(A)
    
    #most_frequent fill benchmark
    def most_frequent_fill(A):
        imp =Imputer(missing_values='NaN', strategy='most_frequent', axis=0, verbose=0, copy=True)
        imp.fit(A)
        return imp.transform(A)
    
    def xrange(x):
        return iter(range(x))
