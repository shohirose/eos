from eos import CubicEosBase, PengRobinsonEos, VanDerWaalsEos, SoaveRedlichKwongEos
from scipy.constants import gas_constant
import pytest


def test_cubic_eos_base():
    pc = 1e6
    tc = 300.0
    omega_a = 1.0
    omega_b = 1.0
    eos = CubicEosBase(pc, tc, omega_a, omega_b)

    assert eos.critical_pressure == pc
    assert eos.critical_temperature == tc

    R = gas_constant
    assert eos._calc_attraction_param() == omega_a*(R*tc)**2/pc
    assert eos._calc_repulsion_param() == omega_b*R*tc/pc

    p = 2e6
    t = 350
    pr = p/pc
    tr = t/tc
    assert eos.reduced_pressure(p) == pr
    assert eos.reduced_temperature(t) == tr
    assert eos._calc_reduced_attraction_param(p, t) == omega_a*pr/tr**2
    assert eos._calc_reduced_repulsion_param(p, t) == omega_b*pr/tr


def test_van_der_waals_eos():
    # Methane
    pc = 4e6       # Critical pressure [Pa]
    tc = 190.6     # Critical temperature [K]
    eos = VanDerWaalsEos(pc, tc)

    p = 3e6
    t = 180
    state = eos.create_state(p, t)
    z = eos.zfactors(state)
    assert z == pytest.approx([0.207498, 0.275339, 0.616434], 1e-4)
    phi = eos.fugacity_coeff(z, state)
    assert phi == pytest.approx([0.756747, 0.758617, 0.741050], 1e-4)


def test_soave_redlich_kwong_eos():
    # Methane
    pc = 4e6       # Critical pressure [Pa]
    tc = 190.6     # Critical temperature [K]
    omega = 0.008  # Acentric factor
    eos = SoaveRedlichKwongEos(pc, tc, omega)

    assert eos.acentric_factor == omega

    p = 3e6
    t = 180
    state = eos.create_state(p, t)
    z = eos.zfactors(state)
    assert z == pytest.approx([0.152443, 0.310673, 0.536884], 1e-4)
    phi = eos.fugacity_coeff(z, state)
    assert phi == pytest.approx([0.69289, 0.70862, 0.70353], 1e-4)


def test_peng_robinson_eos():
    # Methane
    pc = 4e6       # Critical pressure [Pa]
    tc = 190.6     # Critical temperature [K]
    omega = 0.008  # Acentric factor
    eos = PengRobinsonEos(pc, tc, omega)

    assert eos.acentric_factor == omega

    p = 3e6
    t = 180
    state = eos.create_state(p, t)
    z = eos.zfactors(state)
    assert z == pytest.approx([0.135628, 0.292355, 0.510231], 1e-4)
    phi = eos.fugacity_coeff(z, state)
    assert phi == pytest.approx([0.67210, 0.68819, 0.68362], 1e-4)
