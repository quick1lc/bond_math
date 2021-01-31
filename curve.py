# TODO: add ability to pick discounting method (annual, semi_annual, etc.)
# TODO: code comments
# TODO: add ability to get horizon_grid_dict
# TODO: add the ability to set the order of magnitude for the rates

import numpy as np
import pandas as pd

class curve():
    "Curve Class"

    def __init__(self, spot_dict='EMPTY', term_vector='EMPTY',
                 spot_vector='EMPTY', spot_series='EMPTY'):
        """
        The init
        """

        # Absorb inputs; create or parse spot rate series
        if spot_dict != 'EMPTY':
            self.input_terms = [float(str(t)) for t in spot_dict]
            self.input_spots = [float(str(spot_dict[t])) for t in spot_dict]
            self.spot_series = pd.Series(self.input_spots, index=self.input_terms)

        elif (term_vector != 'EMPTY') and (spot_vector != 'EMPTY'):
            self.input_terms = [float(t) for t in term_vector]
            self.input_spots = [float(s) for s in spot_vector]
            self.spot_series = pd.Series(self.input_spots, index=self.input_terms)

        elif spot_series != 'EMPTY':
            self.input_terms = [float(t) for t in list(spot_series.index)]
            self.input_spots = [float(s) for s in list(spot_series)]

        else:
            raise ValueError('Please provide a spot dictonary (spot_dict), '
                             'a spot series (spot_series), '
                             'or term and spot vectors (term_vector and spot_vector)')

        self.spot_length = len(self.spot_series)
        self.terms = self.input_terms
        self.spots = self.input_spots
        self.discount_factors = 'EMPTY'

    # Helper Functions
    def _frange(self, start, stop, step):
        """
        Create a range with float values and float step size
        """

        # Handle basis points (hundredths of a percent)
        oom = 10000
        start = int(start * oom)
        stop = int(stop * oom)
        step = int(step * oom)

        return [i/oom for i in range(start, stop, step)]

    def _calc_diccount_factors(self, spot_series):
        """
        Calculate the discount factors for a given spot curve

        Currently assumes semi-annual
        """

        discounts = pd.Series()
        discounts = (1+(spot_series/200))**((-2*spot_series.index)/12)
        discounts.at[0] = np.nan

        return discounts

    def _individual_forward_rate(self, future_disc, prev_disc, forward_term):
        """
        Calculate an individual implied forward rate
        """

        return 200*((future_disc/prev_disc)**(1/(2*(-forward_term)/12))-1)

    def fill_curve(self, spot_min_term=None, spot_min_val=None,
                   spot_max_term=None, spot_max_val=None):
        """
        Expand a given spot curve to include individual maturities

        Uses linear interpolation
        """

        # Save original spot curve
        self.orig_spot_series = self.spot_series

        # Set Expanded Start and End
        if str(spot_min_term) == 'None':
            spot_min_term = min(self.spot_series.index)
        if str(spot_min_val) == 'None':
            spot_min_val = min(self.spot_series)
        if str(spot_max_term) == 'None':
            spot_max_term = max(self.spot_series.index)
        if str(spot_max_val) == 'None':
            spot_max_val = max(self.spot_series)

        # Expand Spot Series
        if spot_min_term < min(self.spot_series.index):
            self.spot_series.at[spot_min_term] = spot_min_val
        if spot_max_term < max(self.spot_series.index):
            self.spot_series.at[spot_max_term] = spot_max_val
        self.spot_series = self.spot_series.reindex(self._frange(spot_min_term, spot_max_term+1, 1)).interpolate(method='index')

        # Update Object after Expanding Spot Series
        self.terms = list(self.spot_series.index)
        self.spots = list(self.spot_series)
        self.spot_length = len(self.spot_series)

    def add_discount_factors(self):
        """
        Add discount factors to the curve object

        This is a seperate method (not immediatly called by the init) because
        the class includes a simple method for filling the curve if only specific
        spots are given at instantiation. As such, the use must actively add
        the discount factors after they are sure they have a full curve

        Currently assumes semi-annual
        """

        # Add discount factors to object
        self.discount_factors = self._calc_discount_factors(spot_series=self.spot_series)

    def calc_forward_rates(self, forward_term, numb_of_horizons):
        """
        Calculate implied forward rates for a given term, given a spot curve
        """

        # Check that discount factors have been calculated; error if not
        if self.discount_factors == 'EMPTY':
            raise TypeError('Internal discount_factors is default string EMPTY; '
                            'cannot calculate a horizon grid without first '
                            'adding discount factors')
                            
        # Add spot curve to forwards as horizon 0
        fowrard_series = pd.Series()
        forward_series.at[0] = self.spot_series.at[float(forward_term)]

        # Calculate remaining forward horizons
        for i in range(1, numb_of_horizons+1):
            idx=float(i)

            # Test that you have enough spots to calc a forward rate for this horizon
            if (i+forward_term) > self.spot_length:
                curr_forward = np.nan
            else:
                prev_disc = self.factor_series.at[idx]
                future_disc = self.factor_series.at[idx+forward_term]
                curr_forward = self._individual_forward_rate(future_disc=future_disc,
                                                             prev_disc=prev_disc,
                                                             forward_term=forward_term)

            # Add forward to output
            forward_series.at[idx] = curr_forward

        return forward_series

    def calc_horizon_curve(self, horizon_month):
        """
        Calculate an implied spot curve for a future horizon
        """

        # Check that discount factors have been calculated; error if not
        if self.discount_factors == 'EMPTY':
            raise TypeError('Internal discount_factors is default string EMPTY; '
                            'cannot calculate a horizon grid without first '
                            'adding discount factors')

        # Set up data holders
        horizon_term = []
        horizon_spot = []

        # Loop through terms
        start_idx = self.terms.index(horizon_month) + 1 # Add one because discounts has a None at start
        for term in self.terms:

            future_term = horizon_month + term
            if future_term > max(self.terms):
                continue

            future_idx = self.term.index(future_term) + 1 # Add one because discounts has a None at start
            start_disc = self.discount_factors[start_idx]
            future_disc = self.discount_factors[future_idx]
            curr_spot = self._individual_forward_rate(future_disc=future_disc,
                                                      prev_disc=prev_disc,
                                                      forward_term=forward_term)
            horizon_term.append(term)
            horizon_spot.append(curr_spot)

        # Build horizon spot curve as a series
        output_curve = pd.Series(self.horizon_spot, index=self.horizon_term)

        return output_curve

    def calc_horizon_grid(self, numb_of_horizons, max_term=None):
        """
        Calculate a full implied horizon grid
        """

        # Check that discount factors have been calculated; error if not
        if self.discount_factors == 'EMPTY':
            raise TypeError('Internal discount_factors is default string EMPTY; '
                            'cannot calculate a horizon grid without first '
                            'adding discount factors')

        # Get max term if not supplied
        if str(max_term) == 'None':
            max_term = max(self.spot_series.index)

        # Build Grid
        horizon_grid_dict = {}
        temp_df = pd.Dataframe()
        for term in list(self.spot_series.index):

            if t > max_term:
                continue

            self.horizon_grid[term] = self.calc_forward_rates(forward_term=term,
                                                              numb_of_horizons=numb_of_horizons)
            temp_df[term] = self.horizon_grid[term]

        # Transpose DataFrame
        horizon_grid_df = temp_df.T

        return horizon_grid_df
