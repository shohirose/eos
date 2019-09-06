import numpy as np
from abc import ABCMeta, abstractmethod
import math
from scipy.constants import gas_constant


class CubicEOS(metaclass=ABCMeta):
    """
    Base class for cubic equation of states.
    """

    def __init__(self, pc, tc, omega, a, b):
        """
        Parameters
        ----------
        pc : float
            Critical pressure
        tc : float
            Critical temperature
        omega : float
            Accentric factor
        a : float
            Attraction parameter without temperature correction
        b : float
            Volume parameter
        """
        self.__pc = pc
        self.__tc = tc
        self.__omega = omega
        self.__a = a
        self.__b = b

    @property
    def pc(self):
        """ Critical pressure """
        return self.__pc

    @property
    def tc(self):
        """ Critical temperature """
        return self.__tc
    
    @property
    def omega(self):
        """ Accentric factor """
        return self.__omega
    
    @property
    def a(self):
        """ Attraction parameter without temperature correction """
        return self.__a

    @property
    def b(self):
        """ Volume parameter """
        return self.__b

    @staticmethod
    def _calc_a(pc, tc, f):
        """
        Computes the attraction parameter without temperature correction

        Parameters
        ----------
        pc : float
            Critical pressure
        tc : float
            Critical temperature
        f : float
            Coefficient

        Returns
        -------
        float
            Attraction parameter
        """
        return f*(gas_constant*tc)**2/pc

    @staticmethod
    def _calc_b(pc, tc, f):
        """
        Computes the volume parameter

        Parameters
        ----------
        pc : float
            Critical pressure
        tc : float
            Critical temperature
        f : float
            Coefficient
            
        Returns
        -------
        float
            Volume parameter
        """
        return f*gas_constant*tc/pc

    @staticmethod
    def _calc_ar(pr, tr, alpha, f):
        """
        Computes A = a*P/(R*T)^2

        Parameters
        ----------
        pr : float
            Reduced pressure
        tr : float
            Reduced temperature
        alpha : float
            Temperature correction factor
        f : float
            Coefficient

        Returns
        -------
        float
            A = a*P/(R*T)^2
        """
        return f*alpha*pr/tr**2

    @staticmethod
    def _calc_br(pr, tr, f):
        """
        Computes B = b*P/(R*T)

        Parameters
        ----------
        pr : float
            Reduced pressure
        tr : float
            Reduced temperature
        f : float
            Coefficient

        Returns
        -------
        float
            B = b*P/(R*T)
        """
        return f*pr/tr

    @abstractmethod
    def calc_pressure(self, t, v):
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
        pass

    @abstractmethod
    def set(self, p, t):
        """
        Computes parameters required for calc_zfactor() and
        calc_fugacity_coeff().

        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature
        """
        pass

    @abstractmethod
    def calc_zfactor(self):
        """
        Computes Z-factor.
        set() must be called in advance.

        Returns
        -------
        numpy.array
            Z-factors
        """
        pass

    @abstractmethod
    def calc_fugacity_coeff(self, z):
        """
        Computes fugacity coefficient.
        set() must be called in advance.

        Parameters
        ----------
        z : float
            Z-factor

        Returns
        -------
        float
            Fugacity coefficient
        """
        pass


class VanDerWaalsEOS(CubicEOS):
    """
    van der Waals equation of state.
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
        a = self._calc_a(pc, tc, self.__OMEGA_A)
        b = self._calc_b(pc, tc, self.__OMEGA_B)
        super().__init__(pc, tc, None, a, b)

    def calc_pressure(self, t, v):
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
        return gas_constant*t/(v - self.b) - self.a/v**2

    def set(self, p, t):
        """
        Computes parameters required for calc_zfactor() and
        calc_fugacity_coeff().

        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature
        """
        pr = p/self.pc
        tr = t/self.tc
        self.__ar = self._calc_ar(pr, tr, 1.0, self.__OMEGA_A)
        self.__br = self._calc_br(pr, tr, self.__OMEGA_B)

    def calc_zfactor(self):
        """
        Computes Z-factor.
        set() must be called in advance.

        Returns
        -------
        numpy.array
            Z-factors
        """
        a = self.__ar
        b = self.__br
        p = [1.0, -b - 1.0, a, -a*b]
        x = np.roots(p)
        # Returns only real roots
        return np.array([xi.real for xi in x if not np.iscomplex(xi)])

    def calc_fugacity_coeff(self, z):
        """
        Computes fugacity coefficient.
        set() must be called in advance.

        Parameters
        ----------
        z : float
            Z-factor

        Returns
        -------
        float
            Fugacity coefficient
        """
        a = self.__ar
        b = self.__br
        return np.exp(-np.log(z - b) - a/z + z - 1.0)


class SoaveRedlichKwongEOS(CubicEOS):
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
            Accentric factor
        """
        a = self._calc_a(pc, tc, self.__OMEGA_A)
        b = self._calc_b(pc, tc, self.__OMEGA_B)
        super().__init__(pc, tc, omega, a, b)
        self.__m = self.__calc_m(omega)

    @property
    def m(self):
        """ Accentric correction factor """
        return self.__m

    @staticmethod
    def __calc_m(omega):
        """
        Computes the correction factor due to accentric factor.

        Parameters
        ----------
        omega : float
            Accentric factor

        Returns
        -------
        float
            The correction factor due to accentric factor.
        """
        return 0.48 + 1.574*omega - 0.176*omega**2

    @staticmethod
    def __calc_alpha(tr, m):
        """
        Computes temperature correction factor.

        Parameters
        -----------
        tr : float
            Reduced temperature
        m : float
            Accentric correction factor

        Returns
        -------
        float
            Temperature correction factor
        """
        return (1.0 + m*(1.0 - np.sqrt(tr)))**2

    def calc_pressure(self, t, v):
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
        a = self.a
        b = self.b
        tr = t/self.tc
        alpha = self.__calc_alpha(tr, self.m)
        return gas_constant*t/(v - b) - a*alpha/(v*(v + b))

    def set(self, p, t):
        """
        Computes parameters required for calc_zfactor() and
        calc_fugacity_coeff().

        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature
        """
        pr = p/self.pc
        tr = t/self.tc
        alpha = self.__calc_alpha(tr, self.m)
        self.__ar = self._calc_ar(pr, tr, alpha, self.__OMEGA_A)
        self.__br = self._calc_br(pr, tr, self.__OMEGA_B)

    def calc_zfactor(self):
        """
        Computes Z-factor.
        set() must be called in advance.

        Returns
        -------
        numpy.array
            Z-factors
        """
        a = self.__ar
        b = self.__br
        p = [1.0, -1.0, a - b*(1.0 + b), -a*b]
        x = np.roots(p)
        # Returns only real roots
        return np.array([xi.real for xi in x if not np.iscomplex(xi)])

    def calc_fugacity_coeff(self, z):
        """
        Computes fugacity coefficient.
        set() must be called in advance.

        Parameters
        ----------
        z : float
            Z-factor

        Returns
        -------
        float
            Fugacity coefficient
        """
        a = self.__ar
        b = self.__br
        return np.exp(z - 1.0 - np.log(z - b) - a**2/b*np.log(b/z + 1.0))


class PengRobinsonEOS(CubicEOS):
    """
    Peng-Robinson equation of state.
    """
    __SQRT2 = math.sqrt(2.0)
    __DELTA1 = 1.0 + math.sqrt(2.0)
    __DELTA2 = 1.0 - math.sqrt(2.0)
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
            Accentric factor
        """
        a = self._calc_a(pc, tc, self.__OMEGA_A)
        b = self._calc_b(pc, tc, self.__OMEGA_B)
        super().__init__(pc, tc, omega, a, b)
        self.__m = self.__calc_m(omega)

    @staticmethod
    def __calc_m(omega):
        """
        Computes the correction factor due to accentric factor.

        Parameters
        ----------
        omega : float
            Accentric factor

        Returns
        -------
        float
            The correction factor due to accentric factor.
        """
        return 0.3796 + 1.485*omega - 0.1644*omega**2 + 0.01667*omega**3

    @staticmethod
    def __calc_alpha(tr, m):
        """
        Computes temperature correction factor.

        Parameters
        -----------
        tr : float
            Reduced temperature
        m : float
            Accentric correction factor

        Returns
        -------
        float
            Temperature correction factor
        """
        return (1.0 + m*(1.0 - np.sqrt(tr)))**2

    def calc_pressure(self, t, v):
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
        a = self.a
        b = self.b
        tr = t/self.tc
        alpha = self.__calc_alpha(tr, self.__m)
        return gas_constant*t/(v - b) - a*alpha/((v - b)*(v + b) + 2.0*b*v)

    def set(self, p, t):
        """
        Computes parameters required for calc_zfactor() and
        calc_fugacity_coeff().

        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature
        """
        pr = p/self.pc
        tr = t/self.tc
        alpha = self.__calc_alpha(tr, self.__m)
        self.__ar = self._calc_ar(pr, tr, alpha, self.__OMEGA_A)
        self.__br = self._calc_br(pr, tr, self.__OMEGA_B)

    def calc_zfactor(self):
        """
        Computes Z-factor.
        set() must be called in advance.

        Returns
        -------
        numpy.array
            Z-factors
        """
        a = self.__ar
        b = self.__br
        p = [1.0, b - 1.0, a - (3.0*b + 2.0)*b, (-a + (1.0 + b)*b)*b]
        x = np.roots(p)
        # Returns only real roots
        return np.array([xi.real for xi in x if not np.iscomplex(xi)])

    def calc_fugacity_coeff(self, z):
        """
        Computes fugacity coefficient.
        set() must be called in advance.

        Parameters
        ----------
        z : float
            Z-factor

        Returns
        -------
        float
            Fugacity coefficient
        """
        a = self.__ar
        b = self.__br
        SQRT2 = self.__SQRT2
        DELTA1 = self.__DELTA1
        DELTA2 = self.__DELTA2
        return np.exp(z - 1.0 - np.log(z - b) - a/(2.0*SQRT2*b)
                      * np.log((z + DELTA1*b)/(z - DELTA2*b)))


def create(kind, pc, tc, omega):
    """
    Create an EoS object.

    Parameters
    ----------
    kind : string
        Type of EoS
    pc : float
        Critical pressure
    tc : float
        Critical temperature
    omega : float
        Accentric factor

    Returns
    -------
    VanDerWaalsEOS, SoaveRedlichKwongEOS, or PengRobinsonEOS
    """
    if kind == 'VDW':
        return VanDerWaalsEOS(pc, tc)
    elif kind == 'SRK':
        return SoaveRedlichKwongEOS(pc, tc, omega)
    elif kind == 'PR':
        return PengRobinsonEOS(pc, tc, omega)
    else:
        raise ValueError('kind must be VDW, SRK, or PR.')
