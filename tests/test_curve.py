import unittest
import pandas as pd
from bond_math import curve

term = [1,3,5,7,10]
spot = [1,2,3,4,5]
# series = pd.Series(data=spot, index=term)

class curve_test(unittest.TestCase):
    def test_frange(self):
        c = curve(term_vector=term, spot_vector=spot)
        self.assertEqual(c._frange(1.5,4.0,0.5),[1.5,2.0,2.5,3.0,3.5])
        self.assertEqual(c._frange(2,4,1),[2,3])
        self.assertNotEqual(c._frange(1.5,4.0,0.5),[2,3,4])
        self.assertNotEqual(c._frange(2,4,1),[2,3,4])

    # def test_calc_discount_factors(self):
    #     c = curve(term_vector=term, spot_vector=spot)
    #     pd.testing.assert_series_equal(c._calc_discount_factors(series, 2, 100),series)

    # def test_add_discount_factors(self):
    #     c = curve(term_vector=term, spot_vector=spot)
