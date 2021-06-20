function vm = molarVolume(P,T,z)
% Calculate molar volume
%
% Parameters
% ----------
% P : Pressure [Pa]
% T : Temperature [K]
% z : Z-factor
%
% Returns
% -------
% vm : Molar volume [m3/mol]
R = eos.ThermodynamicConstants.Gas;
vm = z*R.*T./P;
end