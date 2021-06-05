
from scipy.constants import gas_constant
import numpy as np

class CubicEosBase:
    """
    Base class for cubic equation of states.
    """

    def __init__(self, pc, tc, omega_a, omega_b):
        """
        Parameters
        ----------
        pc : float
            Critical pressure
        tc : float
            Critical temperature
        omega_a : float
            Coefficient for attraction parameter
        omega_b : float
            Coefficient for repulsion parameter
        """
        self.__pc = pc
        self.__tc = tc
        self.__omega_a = omega_a
        self.__omega_b = omega_b
        self.__a = self._calc_attraction_param()
        self.__b = self._calc_repulsion_param()

    @property
    def critical_pressure(self):
        return self.__pc

    @property
    def critical_temperature(self):
        return self.__tc

    @property
    def attraction_param(self):
        return self.__a

    @property
    def repulsion_param(self):
        return self.__b

    def reduced_pressure(self, p):
        """
        Computes reduced pressure

        Parameters
        ----------
        p : float
            Pressure

        Returns
        -------
        float
            Reduced pressure
        """
        return p/self.__pc

    def reduced_temperature(self, t):
        """
        Computes reduced temperature

        Parameters
        ----------
        t : float
            Temperature

        Returns
        -------
        float
            Reduced temperature
        """
        return t/self.__tc

    def _calc_attraction_param(self):
        """ Computes attraction parameter """
        return self.__omega_a*(gas_constant*self.__tc)**2/self.__pc

    def _calc_repulsion_param(self):
        """ Computes repulsion parameter """
        return self.__omega_b*gas_constant*self.__tc/self.__pc

    def _calc_reduced_attraction_param(self, p, t, alpha=1.0):
        """
        Computes A = alpha*a*P/(R*T)^2

        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature
        alpha : float
            Temperature correction factor. The default value is 1.0.

        Returns
        -------
        float
            A = alpha*a*P/(R*T)^2
        """
        pr = p/self.__pc
        tr = t/self.__tc
        return self.__omega_a*alpha*pr/tr**2

    def _calc_reduced_repulsion_param(self, p, t):
        """
        Computes B = b*P/(R*T)

        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature

        Returns
        -------
        float
            B = b*P/(R*T)
        """
        pr = p/self.__pc
        tr = t/self.__tc
        return self.__omega_b*pr/tr

    @staticmethod
    def solve_cubic_eq(eq):
        """
        Solves a cubic equation

        Parameters
        ----------
        eq : list
            A list of coefficients of a cubic equation
        
        Returns
        -------
        numpy.array
            Real roots of a cubic equation in the ascending order
        """
        x = np.roots(eq)
        return np.sort(np.array([xi.real for xi in x if not np.iscomplex(xi)]))