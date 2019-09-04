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
        self.pc = pc
        self.tc = tc

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

    Parameters
    ----------
    pc : double
        Critical pressure
    tc : double
        Critical temperature
    a : double
        Attraction parameter
    b : double
        Volume parameter
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
        super().__init__(pc, tc)
        self.a = self.calc_a(pc, tc)
        self.b = self.calc_b(pc, tc)

    @staticmethod
    def calc_a(pc, tc):
        """
        Computes the attraction parameter

        Parameters
        ----------
        pc : double
            Critical pressure
        tc : double
            Critical temperature

        Returns
        -------
        double
            Attraction parameter
        """
        return 0.421875*(gas_constant*tc)**2/pc

    @staticmethod
    def calc_b(pc, tc):
        """
        Computes the volume parameter

        Parameters
        ----------
        pc : double
            Critical pressure
        tc : double
            Critical temperature

        Returns
        -------
        double
            Volume parameter
        """
        return 0.125*gas_constant*tc/pc

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

    @staticmethod
    def __calc_A(pr, tr):
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
        return 0.421875*pr/tr**2

    @staticmethod
    def __calc_B(pr, tr):
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
        return 0.125*pr/tr

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
        self.__A = self.__calc_A(pr, tr)
        self.__B = self.__calc_B(pr, tr)

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
    pc : double
        Critical pressure
    tc : double
        Critical temperature
    omega : double
        Accentric factor
    a : double
        Attraction parameter
    b : double
        Volume parameter
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
        super().__init__(pc, tc)
        self.omega = omega
        self.a = self.calc_a(pc, tc)
        self.b = self.calc_b(pc, tc)

    @staticmethod
    def calc_a(pc, tc):
        """
        Computes attraction parameter

        Returns
        -------
        double
            Attraction parameter
        """
        return 0.42748*(gas_constant*self.tc)**2/self.pc

    @staticmethod
    def calc_b(pc, tc):
        """
        Computes volume parameter

        Returns
        -------
        double
            Volume parameter
        """
        return 0.08664*gas_constant*self.tc/self.pc

    @staticmethod
    def calc_alpha(omega, tr):
        """
        Computes temperature correction factor.

        Parameters
        -----------
        omega : double
            Accentric factor
        tr : double
            Reduced temperature

        Returns
        -------
        double
            Temperature correction factor
        """
        m = 0.48 + 1.574*omega - 0.176*omega**2
        return (1.0 + m*(1.0 - np.sqrt(tr)))**2

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
        alpha = self.calc_alpha(self.omega, tr)
        return gas_constant*t/(v - b) - a*alpha/(v*(v + b))

    @staticmethod
    def __calc_A(pr, tr, alpha):
        """
        Computes A = a*P/(R*T)^2

        Parameters
        ----------
        pr : double
            Reduced pressure
        tr : double
            Reduced temperature
        alpha : double
            Temperature correction factor

        Returns
        -------
        double
            A = a*P/(R*T)^2
        """
        return 0.42748*alpha*pr/tr**2

    @staticmethod
    def __calc_B(pr, tr):
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
        return 0.08664*pr/tr

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
        alpha = self.calc_alpha(self.omega, tr)
        self.__A = self.__calc_A(pr, tr, alpha)
        self.__B = self.__calc_B(pr, tr)

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
    pc : double
        Critical pressure
    tc : double
        Critical temperature
    omega : double
        Accentric factor
    a : double
        Attraction parameter
    b : double
        Volume parameter
    """
    SQRT2 = math.sqrt(2.0)

    def __init__(self, pc, tc, omega):
        super().__init__(pc, tc)
        self.omega = omega
        self.a = self.calc_a(pc, tc)
        self.b = self.calc_b(pc, tc)

    @staticmethod
    def calc_a(pc, tc):
        """
        Computes attraction parameter

        Parameters
        ----------
        pc : double
            Critical pressure
        tc : double
            Critical temperature

        Returns
        ------
        double
            Attraction parameter
        """
        return 0.45724*(gas_constant*tc)**2/pc

    @staticmethod
    def calc_b(pc, tc):
        """
        Computes volume parameter

        Parameters
        ----------
        pc : double
            Critical pressure
        tc : double
            Critical temperature

        Returns
        ------
        double
            Volume parameter
        """
        return 0.07780*gas_constant*tc/pc

    @staticmethod
    def calc_alpha(omega, tr):
        """
        Computes temperature correction factor.

        Parameters
        -----------
        omega : double
            Accentric factor
        tr : double
            Reduced temperature

        Returns
        -------
        double
            Temperature correction factor
        """
        m = 0.3796 + 1.485*omega - 0.1644*omega**2 + 0.01667*omega**3
        return (1.0 + m*(1.0 - np.sqrt(tr)))**2

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
        alpha = self.calc_alpha(self.omega, tr)
        return gas_constant*t/(v - b) - a*alpha/((v - b)*(v + b) + 2.0*b*v)

    @staticmethod
    def __calc_A(pr, tr, alpha):
        """
        Computes A = a*P/(R*T)^2

        Parameters
        ----------
        pr : double
            Reduced pressure
        tr : double
            Reduced temperature
        alpha : double
            Temperature correction factor
        """
        return 0.45724*alpha*pr/tr**2

    @staticmethod
    def __calc_B(pr, tr):
        """
        Computes B = b*P/(R*T)

        Parameters
        ----------
        pr : double
            Reduced pressure
        tr : double
            Reduced temperature
        """
        return 0.07780*pr/tr

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
        alpha = self.calc_alpha(self.omega, tr)
        self.__A = self.__calc_A(pr, tr, alpha)
        self.__B = self.__calc_B(pr, tr)

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
        return np.exp(z - 1.0 - np.log(z - B) - A/(2.0*self.SQRT2*B)
                      * np.log((z + (1.0 + self.SQRT2)*B) /
                               (z - (1.0 - self.SQRT2)*B)))
