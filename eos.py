import numpy as np
from abc import ABCMeta, abstractmethod
import math
from scipy.constants import gas_constant


class CubicEOS(metaclass=ABCMeta):
    """
    Base class for cubic equation of states.

    Parameters
    ----------
    pc : double
        Critical pressure
    tc : double
        Critical temperature
    a : double
        Attraction parameter without temperature correction
    b : double
        Volume parameter
    __fa : double
        Coefficient for attraction parameter
    __fb : double
        Coefficient for volume parameter
    """

    def __init__(self, pc, tc, fa, fb):
        """
        Parameters
        ----------
        pc : double
            Critical pressure
        tc : double
            Critical temperature
        fa : double
            Coefficient for attraction parameter
        fb : double
            Coefficient for volume parameter
        """
        self.pc = pc
        self.tc = tc
        self.__fa = fa
        self.__fb = fb
        self.a = self.__calc_a()
        self.b = self.__calc_b()

    def __calc_a(self):
        """
        Computes the attraction parameter

        Returns
        -------
        double
            Attraction parameter without temperature correction
        """
        return self.__fa*(gas_constant*self.tc)**2/self.pc

    def __calc_b(self):
        """
        Computes the volume parameter

        Returns
        -------
        double
            Volume parameter
        """
        return self.__fb*gas_constant*self.tc/self.pc

    def calc_A(self, pr, tr):
        """
        Computes A = a*P/(R*T)^2

        Parameters
        ----------
        pr : double
            Reduced pressure
        tr : double
            Reduced temperature

        Returns
        -------
        double
            A = a*P/(R*T)^2
        """
        return self.__fa*pr/tr**2

    def calc_B(self, pr, tr):
        """
        Computes B = b*P/(R*T)

        Parameters
        ----------
        pr : double
            Reduced pressure
        tr : double
            Reduced temperature

        Returns
        -------
        double
            B = b*P/(R*T)
        """
        return self.__fb*pr/tr

    @abstractmethod
    def calc_pressure(self, t, v):
        """
        Computes pressure.

        Parameters
        ----------
        t : double
            Temperature
        v : double
            Volume

        Returns
        -------
        double
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
        p : double
            Pressure
        t : double
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
        z : double
            Z-factor

        Returns
        -------
        double
            Fugacity coefficient
        """
        pass


class VanDerWaalsEOS(CubicEOS):
    """
    van der Waals equation of state.
    """

    def __init__(self, pc, tc):
        """
        Parameters
        ----------
        pc : double
            Critical pressure
        tc : double
            Critical temperature
        """
        super().__init__(pc, tc, 0.421875, 0.125)

    def calc_pressure(self, t, v):
        """
        Computes pressure.

        Parameters
        ----------
        t : double
            Temperature
        v : double
            Volume

        Returns
        -------
        double
            Pressure
        """
        return gas_constant*t/(v - self.b) - self.a/v**2

    def set(self, p, t):
        """
        Computes parameters required for calc_zfactor() and
        calc_fugacity_coeff().

        Parameters
        ----------
        p : double
            Pressure
        t : double
            Temperature
        """
        pr = p/self.pc
        tr = t/self.tc
        self.__A = self.calc_A(pr, tr)
        self.__B = self.calc_B(pr, tr)

    def calc_zfactor(self):
        """
        Computes Z-factor.
        set() must be called in advance.

        Returns
        -------
        numpy.array
            Z-factors
        """
        A = self.__A
        B = self.__B
        p = [1.0, -B - 1.0, A, -A*B]
        x = np.roots(p)
        # Returns only real roots
        return np.array([xi.real for xi in x if not np.iscomplex(xi)])

    def calc_fugacity_coeff(self, z):
        """
        Computes fugacity coefficient.
        set() must be called in advance.

        Parameters
        ----------
        z : double
            Z-factor

        Returns
        -------
        double
            Fugacity coefficient
        """
        A = self.__A
        B = self.__B
        return np.exp(-np.log(z - B) - A/z + z - 1.0)


class SoaveRedlichKwongEOS(CubicEOS):
    """
    Soave-Redlich-Kwong equation of state.

    Parameters
    ----------
    omega : double
        Accentric factor
    """

    def __init__(self, pc, tc, omega):
        """
        Parameters
        ----------
        pc : double
            Critical pressure
        tc : double
            Critical temperature
        omega : double
            Accentric factor
        """
        super().__init__(pc, tc, 0.42748, 0.08664)
        self.omega = omega
        self.__m = self.__calc_m()

    def __calc_m(self):
        """
        Computes the correction factor due to accentric factor.

        Returns
        -------
        double
            The correction factor due to accentric factor.
        """
        x = self.omega
        return 0.48 + 1.574*x - 0.176*x**2

    def __calc_alpha(self, tr):
        """
        Computes temperature correction factor.

        Parameters
        -----------
        tr : double
            Reduced temperature

        Returns
        -------
        double
            Temperature correction factor
        """
        return (1.0 + self.__m*(1.0 - np.sqrt(tr)))**2

    def calc_pressure(self, t, v):
        """
        Computes pressure.

        Parameters
        ----------
        t : double
            Temperature
        v : double
            Volume

        Returns
        -------
        double
            Pressure
        """
        a = self.a
        b = self.b
        tr = t/self.tc
        alpha = self.__calc_alpha(tr)
        return gas_constant*t/(v - b) - a*alpha/(v*(v + b))

    def set(self, p, t):
        """
        Computes parameters required for calc_zfactor() and
        calc_fugacity_coeff().

        Parameters
        ----------
        p : double
            Pressure
        t : double
            Temperature
        """
        pr = p/self.pc
        tr = t/self.tc
        alpha = self.__calc_alpha(tr)
        self.__A = self.calc_A(pr, tr)*alpha
        self.__B = self.calc_B(pr, tr)

    def calc_zfactor(self):
        """
        Computes Z-factor.
        set() must be called in advance.

        Returns
        -------
        numpy.array
            Z-factors
        """
        A = self.__A
        B = self.__B
        p = [1.0, -1.0, A - B*(1.0 + B), -A*B]
        x = np.roots(p)
        # Returns only real roots
        return np.array([xi.real for xi in x if not np.iscomplex(xi)])

    def calc_fugacity_coeff(self, z):
        """
        Computes fugacity coefficient.
        set() must be called in advance.

        Parameters
        ----------
        z : double
            Z-factor

        Returns
        -------
        double
            Fugacity coefficient
        """
        A = self.__A
        B = self.__B
        return np.exp(z - 1.0 - np.log(z - B) - A**2/B*np.log(B/z + 1.0))


class PengRobinsonEOS(CubicEOS):
    """
    Peng-Robinson equation of state.

    Parameters
    ----------
    SQRT2 : double
        the square root of 2.
    omega : double
        Accentric factor
    """
    __SQRT2 = math.sqrt(2.0)

    def __init__(self, pc, tc, omega):
        super().__init__(pc, tc, 0.45724, 0.07780)
        self.omega = omega
        self.__m = self.__calc_m()

    def __calc_m(self):
        """
        Computes the correction factor due to accentric factor.

        Returns
        -------
        double
            The correction factor due to accentric factor.
        """
        x = self.omega
        return 0.3796 + 1.485*x - 0.1644*x**2 + 0.01667*x**3

    def __calc_alpha(self, tr):
        """
        Computes temperature correction factor.

        Parameters
        -----------
        tr : double
            Reduced temperature

        Returns
        -------
        double
            Temperature correction factor
        """
        return (1.0 + self.__m*(1.0 - np.sqrt(tr)))**2

    def calc_pressure(self, t, v):
        """
        Computes pressure.

        Parameters
        ----------
        t : double
            Temperature
        v : double
            Volume

        Returns
        -------
        double
            Pressure
        """
        a = self.a
        b = self.b
        tr = t/self.tc
        alpha = self.__calc_alpha(tr)
        return gas_constant*t/(v - b) - a*alpha/((v - b)*(v + b) + 2.0*b*v)

    def set(self, p, t):
        """
        Computes parameters required for calc_zfactor() and
        calc_fugacity_coeff().

        Parameters
        ----------
        p : double
            Pressure
        t : double
            Temperature
        """
        pr = p/self.pc
        tr = t/self.tc
        alpha = self.__calc_alpha(tr)
        self.__A = self.calc_A(pr, tr)*alpha
        self.__B = self.calc_B(pr, tr)

    def calc_zfactor(self):
        """
        Computes Z-factor.
        set() must be called in advance.

        Returns
        -------
        numpy.array
            Z-factors
        """
        A = self.__A
        B = self.__B
        p = [1.0, B - 1.0, A - (3.0*B + 2.0)*B, (-A + (1.0 + B)*B)*B]
        x = np.roots(p)
        # Returns only real roots
        return np.array([xi.real for xi in x if not np.iscomplex(xi)])

    def calc_fugacity_coeff(self, z):
        """
        Computes fugacity coefficient.
        set() must be called in advance.

        Parameters
        ----------
        z : double
            Z-factor

        Returns
        -------
        double
            Fugacity coefficient
        """
        A = self.__A
        B = self.__B
        SQRT2 = self.__SQRT2
        return np.exp(z - 1.0 - np.log(z - B) - A/(2.0*SQRT2*B)
                      * np.log((z + (1.0 + SQRT2)*B) /
                               (z - (1.0 - SQRT2)*B)))
