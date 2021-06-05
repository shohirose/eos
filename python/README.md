# eos

Python package for cubic equation of state (EoS).

This package provides three types of cubic EoS: van der Waals, Soave-Redlich-Kwong, and Peng-Robinson EoS. These are respectively defined as the following classes:

- VanDerWaalsEOS
- SoaveRedlichKwongEOS
- PengRobinsonEOS

The EoS classes have the following methods:

- `pressure(t, v)`: computes pressure at a given temperature `t` and volume `v`.
- `create_state(p, t)`: creates a state at pressure `p` and temperature `t`.
- `zfactors(state)`: computes Z-factors of a given state.
- `ln_fugacity_coeff(z, state)`: computes natural log of fugacity coefficients.
- `fugacity_coeff(z, state)`: computes fugacity coefficients.

`create_state` function creates a state which is an instance of `IsobaricIsothermalState`. This state can be passed to `zfactors`, `ln_fugacity_coeff`, and `fugacity_coeff`.

You can create an EoS by using `create_eos` function.

# Examples

Let us assume the following code is already run.

```python
from eos import create_eos
import numpy as np
import matplotlib.pyplot as plt


# Methane
pc = 4e6       # critical pressure [Pa]
tc = 190.6     # critical temperature [K]
omega = 0.008  # accentric factor
```

Then, Z-factors and fugacity coefficients at a given pressure and temperature can be calculated by using van der Waals EoS in the following way.

```python
# Van der Waals EoS
vdw_eos = create_eos('VDW', pc, tc)
# Creates a state
state = vdw_eos.create_state(p=3.0e6, t=180.0)
# Computes Z-factors
z = vdw_eos.zfactors(state)
# Computes fugacity coefficients
phi = vdw_eos.fugacity_coeff(z, state)
```

Pressure-volume plot at constant temperature can be created in the following mannter.

```python
# Peng-Robinson EoS
pr_eos = create_eos('PR', pc, tc, omega=omega)
# Volume
v = np.logspace(-4.3, -2, 1000)
# Temperature
t = 0.9*tc
# Pressure
p = pr_eos.pressure(t, v)

# Plot
plt.figure()
plt.plot(v, p)
plt.title('PR EoS T/Tc = 0.9')
plt.xlabel('Volume [m3/mol]')
plt.ylabel('Pressure [Pa]')
plt.xscale('log')
plt.ylim(-0.25e7, 1e7)
plt.show()
```
