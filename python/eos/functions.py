from eos import VanDerWaalsEos, PengRobinsonEos, SoaveRedlichKwongEos

def create_eos(kind, pc, tc, **kwargs):
    """
    Create an EoS object.

    Parameters
    ----------
    kind : string
        Type of EoS: VDW, SRK, or PR
    pc : float
        Critical pressure
    tc : float
        Critical temperature
    omega : float
        Acentric factor. Required for SDK and PR.

    Returns
    -------
    VanDerWaalsEos, SoaveRedlichKwongEos, or PengRobinsonEos
    """
    if kind == 'VDW':
        return VanDerWaalsEos(pc, tc)
    elif kind == 'SRK':
        if 'omega' not in kwargs:
            raise KeyError("Key 'omega' is not found!")
        omega = kwargs['omega']
        return SoaveRedlichKwongEos(pc, tc, omega)
    elif kind == 'PR':
        if 'omega' not in kwargs:
            raise KeyError("Key 'omega' is not found!")
        omega = kwargs['omega']
        return PengRobinsonEos(pc, tc, omega)
    else:
        raise ValueError('kind must be VDW, SRK, or PR.')