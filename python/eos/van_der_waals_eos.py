from eos import CubicEosBase, IsobaricIsothermalState
from scipy.constants import gas_constant
import numpy as np

class VanDerWaalsEos(CubicEosBase):
    """
    Van der Waals equation of state.
    """
    __OMEGA_A = 0.421875
    __OMEGA_B = 0.125

    def __init__(self, pc, tc):
        """
        Parameters
        ----------
        pc : float
            Critical pressure
        tc : float
            Critical temperature
        """
        super().__init__(pc, tc, self.__OMEGA_A, self.__OMEGA_B)

    @staticmethod
    def _zfactor_cubic_eq(A, B):
        """
        Computes coefficients of the cubic equation of Z-factor.

        Paramters
        ---------
        A : float
            Reduced attraction parameter
        B : float
            Reduced repulsion parameter

        Returns
        -------
        list
            Coefficients of the cubic equation
        """
        return [1.0, -B - 1.0, A, -A*B]

    @staticmethod
    def _ln_fugacity_coeff_impl(z, A, B):
        """
        Computes the natural log of fugacity coefficient

        Parameters
        ----------
        z : float
            Z-factor
        A : float
            Reduced attraction parameter
        B : float
            Reduced repulsion parameter

        Returns
        -------
        float
            Natural log of fugacity coefficient
        """
        return -np.log(z - B) - A/z + z - 1.0

    def pressure(self, t, v):
        """
        Computes pressure at a given volume and temperature.

        Parameters
        ----------
        t : float
            Temperature
        v : float
            Volume

        Returns
        -------
        float
            Pressure
        """
        a = self.attraction_param
        b = self.repulsion_param
        return gas_constant*t/(v - b) - a/v**2

    def create_state(self, p, t):
        """
        Creates constant pressure-temperature state

        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature

        Returns
        -------
        IsobaricIsothermalState
        """
        A = self._calc_reduced_attraction_param(p, t)
        B = self._calc_reduced_repulsion_param(p, t)
        return IsobaricIsothermalState(p, t, A, B)

    def zfactors(self, state):
        """
        Computes Z-factor

        Parameters
        ----------
        state : IsobaricIsothermalState

        Returns
        -------
        numpy.array
            Z-factors in the ascending order
        """
        A = state.reduced_attraction_param
        B = state.reduced_repulsion_param
        return self.solve_cubic_eq(self._zfactor_cubic_eq(A, B))

    def ln_fugacity_coeff(self, z, state):
        """
        Computes the natural log of fugacity coefficient

        Parameters
        ----------
        z : array_like
            Z-factors
        state : IsobaricIsothermalState

        Returns
        -------
        float
            Natural log of fugacity coefficient
        """
        A = state.reduced_attraction_param
        B = state.reduced_repulsion_param
        return self._ln_fugacity_coeff_impl(z, A, B)

    def fugacity_coeff(self, z, state):
        """
        Computes fugacity coefficient

        Parameters
        ----------
        z : array_like
            Z-factors
        state : IsobaricIsothermalState

        Returns
        -------
        float
            Fugacity coefficient
        """
        return np.exp(self.ln_fugacity_coeff(z, state))
