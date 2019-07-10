import unittest
from deicode.base import _BaseImpute


class Test_BaseImpute(unittest.TestCase):
    """Tests base class imports."""
    def test_no_instantiation(self):
        class Foo_boo(_BaseImpute):
            pass
