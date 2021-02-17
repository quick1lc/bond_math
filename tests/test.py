"""
Functional Test for curve class
"""

import bond_math as bmath
import pandas as pd

# Get expectations
in_curve = pd.read_csv('test_curve_in.csv')
exp_curve_H6 = pd.read_csv('test_curve_H6.csv')
exp_H6 = pd.Series(exp_curve_H6['spot'].tolist(),
                   index=exp_curve_H6['term'].tolist())
exp_fwd= pd.read_csv('test_fwd_exp.csv')

# Calc Horizon Spot
test_curve = bmath.curve(term_vector=in_curve['term'].tolist(),
                         spot_vector=in_curve['spot'].tolist())
test_curve.add_discount_factors(compound_periods=2)
calc_H6 = test_curve.calc_horizon_curve(horizon_month=6)
df = pd.DataFrame([exp_H6,calc_H6])
df = df.T
df.columns=['exp','calc']
df.head()

# Calc Forward Rates
df2 = pd.DataFrame()
for term in exp_fwd.columns:
    if term == 'Horizon':
        continue

    df2[f'exp_{term}'] = exp_fwd[term]
    calc_fwd = test_curve.calc_forward_rates(forward_term=float(term),
                                             numb_of_horizons=24)

    df2[f'calc_{term}'] = calc_fwd

df2.head()

# Fill Curve
term = [1,3,5,7,10]
spot = [1,2,3,4,5]
test_curve_2 = bmath.curve(term_vector=term,
                           spot_vector=spot)
test_curve_2.spot_series
test_curve_2.fill_curve()
test_curve_2.spot_series
