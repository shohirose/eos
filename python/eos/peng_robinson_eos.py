from eos import CubicEosBase
import numpy as np
from scipy.constants import gas_constant


class IsobaricIsothermalState:
    """ Constant pressure-temperature state for PR EoS """

    def __init__(self, A, B):
        """
        Parameters
        ----------
        A : float
            Reduced attraction parameter
        B : float
            Reduced repulsion parameter
        """
        self.__A = A
        self.__B = B

    @property
    def reduced_attraction_param(self):
        return self.__A

    @property
    def reduced_repulsion_param(self):
        return self.__B


class PengRobinsonEos(CubicEosBase):
    """
    Peng-Robinson equation of state.
    """
    __SQRT2 = np.sqrt(2.0)
    __DELTA1 = 1.0 + np.sqrt(2.0)
    __DELTA2 = 1.0 - np.sqrt(2.0)
    __OMEGA_A = 0.45724
    __OMEGA_B = 0.07780

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
        return 0.3796 + 1.485*omega - 0.1644*omega**2 + 0.01667*omega**3

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
        return [1.0, B - 1.0, A - 3.0*B**2 - 2.0*B, -A*B + B**2 + B**3]

    @classmethod
    def _calc_Q(cls, z, A, B):
        """
        Computes Q

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
            Q
        """
        return A/(2.0*cls.__SQRT2*B)*np.log((z + cls.__DELTA1*B)/(z + cls.__DELTA2*B))

    @classmethod
    def _ln_fugacity_coeff_impl(cls, z, A, B):
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
        return z - 1.0 - np.log(z - B) - cls._calc_Q(z, A, B)

    def pressure(self, t, v):
        """
        Computes pressure.

        Parameters
        ----------
        t : float
            Temperature
        v : array_like
            Volume

        Returns
        -------
        array_like
            Pressure
        """
        tr = self.reduced_temperature(t)
        alpha = self._calc_alpha(tr)
        a = self.attraction_param
        b = self.repulsion_param
        return gas_constant*t/(v - b) - a*alpha/((v - b)*(v + b) + 2.0*b*v)

    def create_state(self, p, t) -> IsobaricIsothermalState:
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
        return IsobaricIsothermalState(A, B)

    @staticmethod
    def zfactors(state: IsobaricIsothermalState):
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
        return CubicEosBase.solve_cubic_eq(PengRobinsonEos._zfactor_cubic_eq(A, B))

    @classmethod
    def ln_fugacity_coeff(cls, z, state: IsobaricIsothermalState):
        """
        Computes the natural log of fugacity coefficient

        Parameters
        ----------
        z : array_like
            Z-factors
        state : IsobaricIsothermalState

        Returns
        -------
        array_like
            Natural log of fugacity coefficients
        """
        A = state.reduced_attraction_param
        B = state.reduced_repulsion_param
        return cls._ln_fugacity_coeff_impl(z, A, B)

    @classmethod
    def fugacity_coeff(cls, z, state: IsobaricIsothermalState):
        """
        Computes fugacity coefficients

        Parameters
        ----------
        z : array_like
            Z-factors
        state : IsobaricIsothermalState

        Returns
        -------
        array_like
            Fugacity coefficients
        """
        return np.exp(cls.ln_fugacity_coeff(z, state))
