import unittest
import numpy as np
import numpy.testing as npt
import pandas.util.testing as pdt
from deicode.preprocessing import rclr,inverse_rclr
from skbio.stats.composition import closure,clr

class Testpreprocessing(unittest.TestCase):
    def setUp(self):

        self.cdata1=np.array([[2, 2, 6],
                                [4, 4, 2]])
        self.cdata2=[[3, 3, 0]
                    , [0, 4, 2]]
        self.true2=np.array([[ 0.0, 0.0,np.nan],
                        [np.nan,0.34657359, -0.34657359]])
                        
        self.bad1=np.array([1, 2, -1])
        self.bad1
        self._rclr=rclr()
        self._inv=inverse_rclr()
        pass

    def test_rclr(self):

        #test clr works the same if there are no zeros
        cmat=self._rclr.fit_transform(self.cdata1)
        npt.assert_allclose(cmat,clr(self.cdata1.copy()))

        #test a case with zeros :) 
        cmat=self._rclr.fit_transform(self.cdata2)
        npt.assert_allclose(cmat,self.true2)

        with self.assertRaises(ValueError):
            self._rclr.fit_transform(self.bad1)
    
    def test_inverse_rclr(self):

        cmat=self._rclr.fit_transform(self.cdata1)
        npt.assert_allclose(closure(self.cdata1),np.around(self._inv.fit_transform(cmat),1))
        #inverse can not take zero, nan, or inf values (value error)

        pass



