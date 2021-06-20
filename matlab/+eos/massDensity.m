function rho = massDensity(P,T,z,Mw)
% Calculate mass density
%
% Parameters
% ----------
% P : Pressure [Pa]
% T : Temperature [K]
% z : Z-factor
% Mw : Molecular weight [g/mol]
%
% Returns
% -------
% rho : Mass density [kg/m3]
R = eos.ThermodynamicConstants.Gas;
rho = P.*Mw*1e-3./(z*R.*T);
end