"""
Functional Test for curve class
"""

from curve import curve
import pandas as pd

# Get expectations
in_curve = pd.read_csv('test_curve_in.csv')
exp_curve_H6 = pd.read_csv('test_curve_H6.csv')
exp_H6 = pd.Series(exp_curve_H6['spot'].tolist(),
                   index=exp_curve_H6['term'].tolist())
exp_fwd= pd.read_csv('test_fwd_exp.csv')

# Calc Horizon Spot
test_curve = curve(term_vector=in_curve['term'].tolist(),
                   spot_vector=in_curve['spot'].tolist())
test_curve.add_discount_factors(compound_periods=2)
calc_H6 = test_curve.calc_horizon_curve(horizon_month=6)
df = pd.DataFrame([exp_H6,calc_H6])
df = df.T
df.columns=['exp','calc']

# Calc Forward Rates
for term in exp_fwd.columns:
    if term == 'Horizon':
        continue

    calc_fwd = test_curve.calc_forward_rates(forward_term=float(term),
                                             numb_of_horizons=24)
    
    df2 = pd.DataFrame([exp_fwd[term],calc_fwd])
    df2 = df2.T
