import eos.purecomp.VanDerWaalsEos
import eos.purecomp.SoaveRedlichKwongEos
import eos.purecomp.PengRobinsonEos

% Methane
Pc = 4e6;       % Critical pressure [Pa]
Tc = 190.6;     % Critical temperature [K]
omega = 0.008;  % Acentric factor
Mw = 16.0425;   % Molecular weight [g/mol]

P = 3.0e6; % Pressure [Pa]
T = 180;   % Temperature [K]

%% VanDerWaalsEos test:
eos = VanDerWaalsEos(Pc,Tc,Mw);
[z,s] = eos.zFactors(P,T);
z = sort(z);
assert(max(abs((z - [0.207498, 0.275339, 0.616434]')./z)) < 1e-3);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((phi - [0.756747, 0.758617, 0.741050]')./phi)) < 1e-3);

%% SoaveRedlichKwongEos test:
eos = SoaveRedlichKwongEos(Pc,Tc,omega,Mw);
[z,s] = eos.zFactors(P,T);
z = sort(z);
assert(max(abs((z - [0.152443, 0.310673, 0.536884]')./z)) < 1e-3);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((phi - [0.69289, 0.70862, 0.70353]')./phi)) < 1e-3);

%% PengRobinsonEos test:
eos = PengRobinsonEos(Pc,Tc,omega,Mw);
[z,s] = eos.zFactors(P,T);
z = sort(z);
assert(max(abs((z - [0.135628, 0.292355, 0.510231]')./z)) < 1e-3);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((phi - [0.67210, 0.68819, 0.68362]')./phi)) < 1e-3);
