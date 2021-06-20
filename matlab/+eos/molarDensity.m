function rho = molarDensity(P,T,z)
% Calculate molar density
%
% Parameters
% ----------
% P : Pressure [Pa]
% T : Temperature [K]
% z : Z-factor
%
% Returns
% -------
% rho : Molar density [mol/m3]
R = eos.ThermodynamicConstants.Gas;
rho = P./(z*R.*T);
end