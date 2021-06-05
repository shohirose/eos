from eos import CubicEosBase, IsobaricIsothermalState
import numpy as np
from scipy.constants import gas_constant


class SoaveRedlichKwongEos(CubicEosBase):
    """
    Soave-Redlich-Kwong equation of state.
    """
    __OMEGA_A = 0.42748
    __OMEGA_B = 0.08664

    def __init__(self, pc, tc, omega):
        """
        Parameters
        ----------
        pc : float
            Critical pressure
        tc : float
            Critical temperature
        omega : float
            Acentric factor
        """
        super().__init__(pc, tc, self.__OMEGA_A, self.__OMEGA_B)
        self.__omega = omega
        self.__m = self._calc_m(omega)

    @property
    def acentric_factor(self):
        return self.__omega

    @staticmethod
    def _calc_m(omega):
        """
        Computes the correction factor

        Parameters
        ----------
        omega : float
            Acentric factor

        Returns
        -------
        float
            The correction factor
        """
        return 0.48 + 1.574*omega - 0.176*omega**2

    def _calc_alpha(self, tr):
        """
        Computes temperature correction factor.

        Parameters
        -----------
        tr : float
            Reduced temperature

        Returns
        -------
        float
            Temperature correction factor
        """
        return (1.0 + self.__m*(1.0 - np.sqrt(tr)))**2

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
        return [1.0, -1.0, A - B - B**2, -A*B]

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
        return z - 1.0 - np.log(z - B) - A/B*np.log(B/z + 1.0)

    def pressure(self, t, v):
        """
        Computes pressure.

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
        tr = self.reduced_temperature(t)
        alpha = self._calc_alpha(tr)
        a = self.attraction_param
        b = self.repulsion_param
        return gas_constant*t/(v - b) - a*alpha/(v*(v + b))

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
        tr = self.reduced_temperature(t)
        alpha = self._calc_alpha(tr)
        A = self._calc_reduced_attraction_param(p, t, alpha)
        B = self._calc_reduced_repulsion_param(p, t)
        return IsobaricIsothermalState(p, t, A, B)

    def zfactors(self, state):
        """
        Computes Z-factors

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
