import eos.multicomp.VanDerWaalsEos
import eos.multicomp.SoaveRedlichKwongEos
import eos.multicomp.PengRobinsonEos

% CH4, C3H8
Pc = [4.6e6, 4.246e6]';     % Critical pressure [Pa]
Tc = [190.6, 369.8]';       % Critical temperature [K]
omega = [0.008, 0.152]';    % Acentric factor
Mw = [16.0425, 44.1]';      % Molecular weight [g/mol]
K = [0, 0.009;
     0.09, 0]; % Binary interaction parameters

P = 3.0e6; % Pressure [Pa]
T = 180;   % Temperature [K]
x = [0.85, 0.15]'; % Overall composition
V = [0.001, 0.01, 0.1]'; % Volumes [m3]

%% Test 1: VanDerWaalsEos
eos = VanDerWaalsEos(Pc,Tc,Mw,K);
[z,s] = eos.zFactors(P,T,x);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - 0.161471)./z)) < 1e-3);
assert(max(abs((phi - [0.836787, 0.079908]')./phi)) < 1e-3);

%% Test 2: SoaveRedlichKwongEos
eos = SoaveRedlichKwongEos(Pc,Tc,omega,Mw,K);
[z,s] = eos.zFactors(P,T,x);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - 0.128659)./z)) < 1e-3);
assert(max(abs((phi - [0.810235, 0.073144]')./phi)) < 1e-3);

%% Test 3: PengRobinsonEos
eos = PengRobinsonEos(Pc,Tc,omega,Mw,K);
[z,s] = eos.zFactors(P,T,x);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - 0.111994)./z)) < 1e-3);
assert(max(abs((phi - [0.688407, 0.123469]')./phi)) < 1e-3);