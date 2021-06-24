classdef SoaveRedlichKwongEos < eos.multicomp.CubicEosBase
    % Soave-Redlich-Kwong equation of state.
    
    properties (SetAccess = private)
        AcentricFactor % Acentric factor
    end
    methods (Static)
        function coeffs = zFactorCubicEq(A,B)
            % Compute coefficients of Z-factor cubic equation.
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
        function lnPhi = lnFugacityCoeffImpl(z,x,A,B,Aij,Bi)
            % Computes natural log of fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factor
            % x : Composition
            % A : Attraction parameter of the mixture
            % B : Repulsion parameter of the mixture
            % Aij : Combined attraction parameter between i and j
            % components
            % Bi : Repulsion parameter of each component
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            arguments
                z (1,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric}
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
                Aij (:,:) {mustBeNumeric}
                Bi (:,1) {mustBeNumeric}
            end
            lnPhi = Bi/B*(z - 1) - log(z - B) ...
                - A/B*log(B/z + 1)*(2*(Aij*x)/A - Bi/B);
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
            R = eos.ThermodynamicConstants.Gas;
            P = R*T./(V - b) - a./(V.*(V + b));
        end
        function coeffs = dPdVPolyEq(T,a,b)
            % Compute the coefficients of the polynomial of dPdV = 0.
            %
            % coeffs = DPDVPOLYEQ(T,a,b)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % a : Attraction parameter
            % b : Repulsion parameter
            %
            % Returns
            % -------
            % coeffs : Coefficients of the polynomial
            R = eos.ThermodynamicConstants.Gas;
            coeffs = [R*T, 2*(b*R*T - a), b^2*R*T + 3*a*b, 0, -a*b^3];
        end
    end
    methods
        function obj = SoaveRedlichKwongEos(Pc,Tc,omega,Mw,K)
            % Construct SRK EoS.
            %
            % obj = SOAVEREDLICHKWONGEOS(Pc,Tc,omega,Mw,K)
            %
            % Parameters
            % ----------
            % Pc    : Critical pressure [Pa]
            % Tc    : Critical temperature [K]
            % omega : Acentric factor
            % Mw    : Molecular weight [g/mol]
            % K     : Binary interaction parameters
            %
            % Returns
            % -------
            % obj : An instance of SRK EoS
            arguments
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric}
            end
            obj@eos.multicomp.CubicEosBase(0.42748,0.08664,Pc,Tc,Mw,K)
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw,K)
            % Set parameters.
            %
            % obj = obj.SETPARAMS(Pc,Tc,omega,Mw,K)
            %
            % Parameters
            % ----------
            % Pc    : Critical pressure [Pa]
            % Tc    : Critical temperature [K]
            % omega : Acentric factor
            % Mw    : Molecular weight [g/mol]
            % K     : Binary interaction parameters
            %
            % Returns
            % -------
            % obj : an instance of SRK EoS
            arguments
                obj {mustBeA(obj,'eos.multicomp.SoaveRedlichKwongEos')}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric}
            end
            obj = setParams@eos.multicomp.CubicEosBase(obj,Pc,Tc,Mw,K);
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
            m = 0.48 + 1.574*omega - 0.176*omega.^2;
            alpha = (1 + m.*(1 - sqrt(Tr))).^2;
        end
    end
end