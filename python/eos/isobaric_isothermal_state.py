class IsobaricIsothermalState:
    """ Constant pressure-temperature state """

    def __init__(self, p, t, A, B):
        """
        Parameters
        ----------
        p : float
            Pressure
        t : float
            Temperature
        A : float
            Reduced attraction parameter
        B : float
            Reduced repulsion parameter
        """
        self.__p = p
        self.__t = t
        self.__A = A
        self.__B = B

    @property
    def pressure(self):
        return self.__p

    @property
    def temperature(self):
        return self.__t

    @property
    def reduced_attraction_param(self):
        return self.__A

    @property
    def reduced_repulsion_param(self):
        return self.__B
