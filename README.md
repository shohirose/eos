# eos
Python codes for cubic equation of state (EoS).

This module provides three types of cubic EoS: van der Waals, Soave-Redlich-Kwong, and Peng-Robinson EoS. These are respectively defined as the following classes:

- VanDerWaalsEOS
- SoaveRedlichKwongEOS
- PengRobinsonEOS

The EoS classes have the following methods:

- `calc_pressure(t, v)`: compute pressure at a given temperature `t` and volume `v`
- `calc_zfactor()`: compute Z-factor at a given pressure `p` and temperature `t`
- `calc_fugacity_coeff(z)`: compute fugacity coefficient for a given Z-factor

Please note that `set(p, t)` must be called before calling`calc_zfactor()` and `calc_fugacity_coeff(z)`.


# Examples

Let us assume the following code is already run.

```python
import eos
import numpy as np
import matplotlib.pyplot as plt


# Methane
# critical pressure [Pa]
pc = 4e6
# critical temperature [K]
tc = 190.6
# accentric factor
omega = 0.008
```

Then, Z-factor and fugacity coefficient at a given pressure and temperature can be calculated by using van der Waals EoS in the following way.

```python
# Van der Waals EoS
vdw_eos = eos.VanDerWaalsEOS(pc, tc)
# set P and T
vdw_eos.set(p=3.0e6, t=180.0)
# Computes Z-factor
z = vdw_eos.calc_zfactor()
# Computes fugacity coefficient
phi = vdw_eos.calc_fugacity_coeff(z[0])
```

P-V plot at constant temperature can be created in the following mannter.

```python
# Peng-Robinson EoS
pr_eos = eos.PengRobinsonEOS(pc, tc, omega)
# Volume
v = np.logspace(-4.3, -2, 1000)
# Temperature
t = 0.9*tc
# Pressure
p = pr_eos.calc_pressure(t, v)

# Plot
plt.figure()
plt.plot(v, p)
plt.title('PR-EoS T/Tc = 0.9')
plt.xlabel('V [m3/mol]')
plt.ylabel('P [Pa]')
plt.xscale('log')
plt.ylim(-0.25e7, 1e7)
plt.show()
```
