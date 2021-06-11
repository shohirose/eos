import eos.VanDerWaalsEos
import eos.SoaveRedlichKwongEos
import eos.PengRobinsonEos

% Methane
Pc = 4e6;       % Critical pressure [Pa]
Tc = 190.6;     % Critical temperature [K]
omega = 0.008;  % Acentric factor

P = 3.0e6; % Pressure [Pa]
T = 180;   % Temperature [K]

%% VanDerWaalsEos test:
vdweos = VanDerWaalsEos(Pc,Tc);
[z,A,B] = vdweos.zFactors(P,T);
z = sort(z);
assert(max(abs((z - [0.207498, 0.275339, 0.616434]')./z)) < 1e-3);
phi = vdweos.fugacityCoeff(z,A,B);
assert(max(abs((phi - [0.756747, 0.758617, 0.741050]')./phi)) < 1e-3);

%% SoaveRedlichKwongEos test:
srkeos = SoaveRedlichKwongEos(Pc,Tc,omega);
[z,A,B] = srkeos.zFactors(P,T);
z = sort(z);
assert(max(abs((z - [0.152443, 0.310673, 0.536884]')./z)) < 1e-3);
phi = srkeos.fugacityCoeff(z,A,B);
assert(max(abs((phi - [0.69289, 0.70862, 0.70353]')./phi)) < 1e-3);

%% PengRobinsonEos test:
preos = PengRobinsonEos(Pc,Tc,omega);
[z,A,B] = preos.zFactors(P,T);
z = sort(z);
assert(max(abs((z - [0.135628, 0.292355, 0.510231]')./z)) < 1e-3);
phi = preos.fugacityCoeff(z,A,B);
assert(max(abs((phi - [0.67210, 0.68819, 0.68362]')./phi)) < 1e-3);
