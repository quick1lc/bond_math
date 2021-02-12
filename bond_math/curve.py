# TODO: add ability to get horizon_grid_dict

import numpy as np
import pandas as pd

class curve():
    """
    Cruve helps with simple operations to 'roll a spot curve' forward
    """

    def __init__(self, spot_dict={}, term_vector=[],
                 spot_vector=[], spot_series=pd.Series,
                 order_of_mag=100):
        """
        The init for the curve class. Can be instantiated with a spot curve
        dictonary, pandas Series, or as rate and term vectors. Once created,
        the curve object will hold the given curve and facilitate the
        calculation of forward rate paths and/or future spot curves.

        :param spot_dict: A spot curve given in the form of a dictonary, where
        the each term is a key, and the corresponding rate is the value
        :type spot_dict: dict
        :param term_vector: The terms from a spot curve as a list (must be
        paired with a spot_vector)
        :type term_vector: list
        :param spot_vector: The spot rates from a spot curve as a list (must be
        paired with a term_vector)
        :type spot_vector: list
        :param spot_series: A spot curve given in the form of a pandas series,
        where the index is the term, and the spot rates are the members of the
        series
        :type spot_series: pd.Series
        :param order_of_mag: The divisor needed to convert the given rates into
        a decimal (where 1% = 0.01). If the rates are given in percentages
        (where 1% = 1), then the order of magnitude needs to be 100
        :type order_of_mag: int
        """

        # Absorb inputs; create or parse spot rate series
        if spot_dict:
            self.input_terms = [float(str(t)) for t in spot_dict]
            self.input_spots = [float(str(spot_dict[t])) for t in spot_dict]
            self.spot_series = pd.Series(self.input_spots, index=self.input_terms)

        elif term_vector and spot_vector:
            self.input_terms = [float(t) for t in term_vector]
            self.input_spots = [float(s) for s in spot_vector]
            self.spot_series = pd.Series(self.input_spots, index=self.input_terms)

        elif not spot_series.empty:
            self.input_terms = [float(t) for t in list(spot_series.index)]
            self.input_spots = [float(s) for s in list(spot_series)]
            self.spot_series = spot_series

        else:
            raise ValueError('Please provide a spot dictonary (spot_dict), '
                             'a spot series (spot_series), '
                             'or term and spot vectors (term_vector and spot_vector)')

        # Set up other holders
        self.spot_length = len(self.spot_series)
        self.terms = self.input_terms
        self.spots = self.input_spots
        self.discount_factors = pd.Series
        self.compound_periods = None
        self.order_of_mag = order_of_mag

    # Helper Functions
    def _frange(self, start, stop, step):
        """
        Create a range with float values and float step size

        :param start: The first value of the desired range
        :type start: float
        :param stop: The stoping point for the range, not included in the
        resulting list
        :type stop: float
        :param step: The desired step between items in the resulting range
        :type step: float

        :returns: A list created from the desired range
        """

        # Handle basis points (hundredths of a percent)
        oom = 10000
        start = int(start * oom)
        stop = int(stop * oom)
        step = int(step * oom)

        return [i/oom for i in range(start, stop, step)]

    def _calc_discount_factors(self, spot_series, compound_periods,
                               order_of_mag):
        """
        Calculate the discount factors for a given spot curve; this method
        assumes that all terms are in months

        :param spot_series: The spot rates to be converted to discount factors;
        the series should have rates as the items and terms as the index
        :type spot_series: pandas Series
        :param compound_periods: The number of compounding periods per year
        :type compound_periods: int or float
        :param order_of_mag: The divisor needed to convert the given rates into
        a decimal (where 1% = 0.01). If the rates are given in percentages
        (where 1% = 1), then the order of magnitude needs to be 100
        :type order_of_mag: int

        :returns: a pandas Series of discount factors
        """
        # Example Calc
        # discounts = (1+(spot_series/200))**((-2*spot_series.index)/12)

        discounts = pd.Series()
        discounts = (1+((spot_series/order_of_mag)/compound_periods))**(-compound_periods*(spot_series.index/12))
        discounts.at[0] = np.nan

        return discounts

    def _individual_forward_rate(self, future_disc, start_disc, forward_term,
                                  compound_periods, order_of_mag):
        """
        Calculate an individual implied forward rate using two discount factors;
        this method assumes that all terms are in months

        :param future_disc: The discount factor for the future horizon; should
        be discount factor for the spot term at the desired horizon + the
        desired forward term. If you are trying to the 6 month forward in month
        4, then the future_disc would be the spot discount factor from term 10
        :type future_disc: float
        :param start_disc: The discount factor for the base horizon; when the
        spot curve has been expanded, the start_disc should be the discount factor
        of the spot term that is equal to the desired horizon
        :type start_disc: float
        :param forward_term: The term of the desired forward rate
        :type forward_term: float
        :param compound_periods: The number of compounding periods per year
        :type compound_periods: int or float
        :param order_of_mag: The divisor needed to convert the given rates into
        a decimal (where 1% = 0.01). If the rates are given in percentages
        (where 1% = 1), then the order of magnitude needs to be 100
        :type order_of_mag: int

        :returns: an individual forward rate as a float
        """
        # Example Calc
        # 200*((future_disc/prev_disc)**(1/(2*(-forward_term)/12))-1)

        return (order_of_mag*compound_periods)*((future_disc/start_disc)**(1/(compound_periods*(-forward_term/12)))-1)


    def  _individual_horizon_spot(self, future_disc, start_disc, term,
                                  compound_periods, order_of_mag):
        """
        Calculate an individual implied horizon spot rate using two discount
        factors; this method assumes that all terms are in months

        :param future_disc: The discount factor for the future horizon; should
        be discount factor for the spot term at the desired horizon + the
        desired spot term. If you are trying to get the 2 month spot rate in
        horizon 6, then the future_disc would be the spot discount factor from
        term 8
        :type future_disc: float
        :param start_disc: The spot discount rate from the desired horizon; If
        you are trying to build the spot curve in horizon 6, then the start_disc
        would be the the spot discount factor from term 6
        :type start_disc: float
        :param term: The spot term being requested
        :type term: float
        :param compound_periods: The number of compounding periods per year
        :type compound_periods: int or float
        :param order_of_mag: The divisor needed to convert the given rates into
        a decimal (where 1% = 0.01). If the rates are given in percentages
        (where 1% = 1), then the order of magnitude needs to be 100
        :type order_of_mag: int

        :returns: an individual forward rate as a float
        """
        # Example Calc
        # 200*((start_disc/future_disc)**(1/(2*(term/12)))-1)

        return (order_of_mag*compound_periods)*((start_disc/future_disc)**(1/(compound_periods*(term/12)))-1)


    def fill_curve(self, spot_min_term=None, spot_min_val=None,
                   spot_max_term=None, spot_max_val=None):
        """
        Expands a given spot curve to include individual maturities; also
        allows for the curve to be expanded to longer or shorter maturities
        than were included in the original spot curve.

        Note:
            This method uses linear interpolation, and is not a suitable
            substitute for a proper model.

        :param spot_min_term: The min desired term of the expanded spot curve;
        if not lower than the min term of the original spot curve, the term will
        be inserted inside the curve
        :type spot_min_term: int or float
        :param spot_min_val: The spot rate at the expanded min term
        :type spot_min_val: int or float
        :param spot_max_term: The max desired term of the expanded spot curve;
        if not bigger than the max term of the original spot curve, the term will
        be inserted inside the curve
        :type spot_max_term: int or float
        :param spot_max_val: The spot rate at the expanded max term
        :type spot_max_val: int or float

        :returns: None
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
        if spot_max_term > max(self.spot_series.index):
            self.spot_series.at[spot_max_term] = spot_max_val
        self.spot_series = self.spot_series.reindex(self._frange(spot_min_term, spot_max_term+1, 1)).interpolate(method='index')

        # Update Object after Expanding Spot Series
        self.terms = list(self.spot_series.index)
        self.spots = list(self.spot_series)
        self.spot_length = len(self.spot_series)

    def add_discount_factors(self, compound_periods=1):
        """
        Add discount factors to the curve object (adds the new series to the
        object as self.discount_factors; nothing is returned)

        This is a seperate method (not immediatly called by the init) because
        the class includes a simple method for filling the curve if only specific
        spots are given at instantiation. As such, the use must actively add
        the discount factors after they are sure they have a full curve

        :param compound_periods: The number of compounding periods per year
        :type compound_periods: int or float

        :returns: None
        """

        # Add discount factors to object
        self.discount_factors = self._calc_discount_factors(spot_series=self.spot_series,
                                                            compound_periods=compound_periods,
                                                            order_of_mag=self.order_of_mag)

        # Capture requested number of compounding periods
        self.compound_periods = compound_periods

    def calc_forward_rates(self, forward_term, numb_of_horizons):
        """
        Calculate implied forward rates for a given term, given a spot curve

        :param forward_term: The desired forwrad term
        :type forward_term: int
        :param numb_of_horizons: The desired length of the forward rate path
        :type numb_of_horizons: int

        :returns: The forward rate path as a pandas Series where the horizon
        is the index and the rates are the values
        """

        # Check that discount factors have been calculated; error if not
        if self.discount_factors.empty:
            raise TypeError('Internal discount_factors are empty; '
                            'cannot calculate a horizon grid without first '
                            'adding discount factors using add_discount_factors')

        # Add spot curve to forwards as horizon 0
        forward_series = pd.Series()
        forward_series.at[0] = self.spot_series.at[float(forward_term)]

        # Calculate remaining forward horizons
        for i in range(1, numb_of_horizons+1):
            idx=float(i)

            # Test that you have enough spots to calc a forward rate for this horizon
            if (i+forward_term) > self.spot_length:
                curr_forward = np.nan
            else:
                start_disc = self.discount_factors.at[idx]
                future_disc = self.discount_factors.at[idx+forward_term]
                curr_forward = self._individual_forward_rate(future_disc=future_disc,
                                                             start_disc=start_disc,
                                                             forward_term=forward_term,
                                                             compound_periods=self.compound_periods,
                                                             order_of_mag=self.order_of_mag)

            # Add forward to output
            forward_series.at[idx] = curr_forward

        return forward_series

    def calc_horizon_curve(self, horizon_month):
        """
        Calculate an implied spot curve for a future horizon

        :param horizon_month: The horizon month for which to calculate an
        implied spot curve
        :type horizon_month: int

        :returns: A horizon spot curve as a pandas series where the term is the
        index and the spot rates are the values
        """

        # Check that discount factors have been calculated; error if not
        if self.discount_factors.empty:
            raise TypeError('Internal discount_factors are empty; '
                            'cannot calculate a horizon grid without first '
                            'adding discount factors using add_discount_factors')

        # Set up data holders
        horizon_term = []
        horizon_spot = []

        # Loop through terms
        start_idx = self.terms.index(horizon_month) + 1 # Add one because discounts has a None at start
        for term in self.terms:

            future_term = horizon_month + term
            if future_term > max(self.terms):
                continue

            future_idx = self.terms.index(future_term) + 1 # Add one because discounts has a None at start
            start_disc = self.discount_factors[start_idx]
            future_disc = self.discount_factors[future_idx]
            curr_spot = self._individual_horizon_spot(future_disc=future_disc,
                                                      start_disc=start_disc,
                                                      term=term,
                                                      compound_periods=self.compound_periods,
                                                      order_of_mag=self.order_of_mag)

            horizon_term.append(term)
            horizon_spot.append(curr_spot)

        # Build horizon spot curve as a series
        output_curve = pd.Series(horizon_spot, index=horizon_term)

        return output_curve

    def calc_horizon_grid(self, numb_of_horizons, max_term=None):
        """
        Calculate a full implied horizon grid.

        :param numb_of_horizons: The desired length of the forward rate paths
        :type numb_of_horizons: int
        :param max_term: The max term to be included in the grid
        :type max_term: int

        :returns: The horizon grid as a pandas Dataframe. The columns of the
        grid are the terms and the index is the horizon month
        """

        # Check that discount factors have been calculated; error if not
        if self.discount_factors.empty:
            raise TypeError('Internal discount_factors are empty; '
                            'cannot calculate a horizon grid without first '
                            'adding discount factors using add_discount_factors')

        # Get max term if not supplied
        if str(max_term) == 'None':
            max_term = max(self.terms)

        # Build Grid
        # horizon_grid_dict = {}
        temp_df = pd.Dataframe()
        for term in self.terms:

            if term > max_term:
                continue

            self.horizon_grid[term] = self.calc_forward_rates(forward_term=term,
                                                              numb_of_horizons=numb_of_horizons)
            temp_df[term] = self.horizon_grid[term]

        # Transpose DataFrame
        horizon_grid_df = temp_df.T

        return horizon_grid_df
