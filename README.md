# eos
Python codes for cubic equation of state.

# Examples

```python
import eos
import numpy as np
import matplotlib.pyplot as plt


# Methane
# critical pressure
pc = 4e6
# critical temperature
tc = 190.6
# accentric factor
omega = 0.008

# ---------------------------------------------------------
# Example of Z-factor and fugacity coefficient calculation
# ---------------------------------------------------------
# Van der Waals EoS
vdw_eos = eos.VanDerWaalsEOS(pc, tc)
# set P and T
eos.set(p=3.0e6, t=180.0)
# Computes Z-factor
z = eos.calc_zfactor()
# Computes fugacity coefficient
phi = eos.calc_fugacity_coeff(z[0])

# ---------------------------------------------------------
# Example of P-V plot at constant temperature
# ---------------------------------------------------------
pr_eos = eos.PengRobinsonEOS(pc, tc, omega)
# Volume
v = np.logspace(-4.3, -2, 1000)
# Temperature
t = 0.9*tc
# Pressure
p = pr_eos.calc_pressure(t, v)

plt.figure()
plt.plot(v, p)
plt.title('PR-EoS T/Tc = 0.9')
plt.xlabel('V [m3/mol]')
plt.ylabel('P [Pa]')
plt.xscale('log')
plt.ylim(-0.25e7, 1e7)
plt.show()
```
