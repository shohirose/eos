classdef SoaveRedlichKwongEos < eos.purecomp.CubicEosBase
    % SoaveRedlichKwongEos Soave-Redlich-Kwong equation of state
    %
    %  This class provides methods to calculate thermodynamic properties
    %  based on Soave-Redlich-Kwong equation of state.
    
    properties (SetAccess = private)
        AcentricFactor % Acentric factor
    end
    methods (Static)
        function coeffs = zFactorCubicEq(A,B)
            % Compute coefficients of Z-factor cubic equation
            %
            % coeffs = ZFACTORCUBICEQ(A,B)
            %
            % Parameters
            % ----------
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            %
            % Returns
            % -------
            % coeffs : Coefficients of the cubic equation of Z-factor
            arguments
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
            end
            coeffs = [1, -1, A - B - B^2, -A*B];
        end
        function coeffs = dPdTPolyEq(T,a,b)
            % Compute coefficients of the polynomial of dPdT = 0.
            %
            % coeffs = DPDTPOLYEQ(T,a,b)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % a : Attraction parameter
            % b : Repulsion parameter
            %
            % Returns
            % -------
            % coeffs : Coefficients of the polynomial of dPdT = 0.
            arguments
                T (1,1) {mustBeNumeric}
                a (1,1) {mustBeNumeric}
                b (1,1) {mustBeNumeric}
            end
            R = eos.ThermodynamicConstants.Gas;
            coeffs = [R*T, 2*(b*R*T - a), b^2*R*T + 3*a*b, 0, -a*b^3];
        end
        function lnPhi = lnFugacityCoeff(z,A,B)
            % Compute natural log of fugacity coefficients
            %
            % lnPhi = LNFUGACITYCOEFF(z,s)
            %
            % Parameters
            % ----------
            % z : Z-factor
            % A : Attraction parameter
            % B : Repulsion parameter
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            arguments
                z (:,1) {mustBeNumeric}
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
            end
            lnPhi = z - 1 - log(z - B) - A/B*log(B./z + 1);
        end
        function P = pressureImpl(T,V,a,b)
            % Compute pressure.
            %
            % P = PRESSUREIMPL(T,V,a,b)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % V : Volume [m3]
            % a : Attraction parameter
            % b : Repulsion parameter
            %
            % Returns
            % -------
            % P : Pressure [Pa]
            arguments
                T {mustBeNumeric}
                V {mustBeNumeric}
                a {mustBeNumeric}
                b {mustBeNumeric}
            end
            R = eos.ThermodynamicConstants.Gas;
            P = R*T./(V - b) - a./(V.*(V + b));
        end
    end
    methods
        function obj = SoaveRedlichKwongEos(Pc,Tc,omega,Mw)
            % Construct SRK EoS
            %
            % obj = SOAVEREDLICHKWONGEOS(Pc,Tc,omega,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            %
            % Returns
            % -------
            % obj : SOAVEREDLICHKWONGEOS
            arguments
                Pc (1,1) {mustBeNumeric}
                Tc (1,1) {mustBeNumeric}
                omega (1,1) {mustBeNumeric}
                Mw (1,1) {mustBeNumeric}
            end
            obj@eos.purecomp.CubicEosBase(0.42748,0.08664,Pc,Tc,Mw)
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw)
            % Set parameters
            %
            % obj = obj.SETPARAMS(Pc,Tc,omega,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            %
            % Returns
            % -------
            % obj : SOAVEREDLICHKWONGEOS
            arguments
                obj {mustBeA(obj,'eos.purecomp.SoaveRedlichKwongEos')}
                Pc (1,1) {mustBeNumeric}
                Tc (1,1) {mustBeNumeric}
                omega (1,1) {mustBeNumeric}
                Mw (1,1) {mustBeNumeric}
            end
            obj = setParams@eos.purecomp.CubicEosBase(obj,Pc,Tc,Mw);
            obj.AcentricFactor = omega;
        end
        function alpha = temperatureCorrectionFactor(obj,Tr)
            % Compute temperature correction factor.
            %
            % alpha = obj.TEMPERATURECORRECTIONFACTOR(Tr)
            %
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % alpha : Temperature correction factor
            omega = obj.AcentricFactor;
            m = 0.48 + 1.574*omega - 0.176*omega^2;
            alpha = (1 + m*(1 - sqrt(Tr))).^2;
        end
    end
end