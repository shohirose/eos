%% Example of VanDerWaalsEos class
import eos.VanDerWaalsEos

%% Create an instance of VDW EOS class
% Methane
Pc = 4e6;       % Critical pressure [Pa]
Tc = 190.6;     % Critical temperature [K]
Mw = 16.0425;   % Molecular weight [g/mol]
vdw = eos.VanDerWaalsEos(Pc,Tc,Mw);

%% Plot isothermal lines
b = vdw.RepulsionParam;
V = logspace(log10(b*1.1), log10(1e-2), 100);
T = [150, 170, Tc, 210];
[T,V] = meshgrid(T,V);
P = vdw.pressure(T,V);

figure;
semilogx(V,P);
legend('T=150','T=170','T=Tc','T=210');
axis([1e-5, 1e-2, -3e6, 1e7]);
title('Van der Waals EoS');
xlabel('Volume [m3]');
ylabel('Pressure [Pa]');

%% Computes Z-factors and fugacity
P = 3e6;
T = 180;
% Computes Z-factors
[z,A,B] = vdw.zFactors(P,T);
% Computes fugacity coefficients
phi = vdw.fugacityCoeff(z,A,B);
% Computes density and volume
rhom = vdw.molarDensity(P,T,z);
rhow = vdw.massDensity(P,T,z);
vm = vdw.molarVolume(P,T,z);