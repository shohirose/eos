classdef PengRobinsonEos < eos.multicomp.CubicEosBase
    % PengRobinsonEos Peng-Robinson equation of state
    %
    %  This class provides methods to calculate thermodynamic properties
    %  based on Peng-Robinson equation of state.
    
    properties (Constant, Access = private)
        Sqrt2 = sqrt(2)
        Delta1 = 1 + sqrt(2)
        Delta2 = 1 - sqrt(2)
    end
    properties (SetAccess = private)
        AcentricFactor % Acentric factor
    end
    methods (Static)
        function coeffs = zFactorCubicEq(A,B)
            % Compute coefficients of the cubic equation of Z-factor
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
            coeffs = [1, B - 1, A - 2*B - 3*B^2, -A*B + B^2 + B^3];
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
            Sqrt2 = eos.multicomp.PengRobinsonEos.Sqrt2;
            Delta1 = eos.multicomp.PengRobinsonEos.Delta1;
            Delta2 = eos.multicomp.PengRobinsonEos.Delta2;
            lnPhi = z - 1 - log(z - B) ...
                - A./(2*Sqrt2*B)*log((z + Delta1*B)./(z + Delta2*B)) ...
                *(2*(Aij*x)/A - Bi/B);
        end
        function P = pressureImpl(T,V,a,b)
            % Compute pressure
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
            P = R*T./(V - b) - a./((V - b).*(V + b) + 2*b*V);
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
            coeffs = [R*T, 4*b*R*T - 2*a, 2*(b^2*R*T + a*b), ...
                2*b^2*(a - 2*b*R*T), b^3*(b*R*T - 2*a)];
        end
    end
    methods
        function obj = PengRobinsonEos(Pc,Tc,omega,Mw,K)
            % Construct Peng-Robinson EoS
            %
            % eos = PENGROBINSONEOS(Pc,Tc,omega,Mw,K)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            % K : Binary interaction parameters
            %
            % Returns
            % -------
            % obj : an instance of PENGROBINSONEOS
            arguments
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric}
            end
            obj@eos.multicomp.CubicEosBase(0.45724,0.07780,Pc,Tc,Mw,K);
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw,K)
            % Set parameters
            %
            % obj = obj.SETPARAMS(Pc,Tc,omega,Mw,K)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            % K : Binary interaction composition
            %
            % Returns
            % -------
            % obj : an instance of PENGROBINSONEOS
            arguments
                obj {mustBeA(obj,'eos.multicomp.PengRobinsonEos')}
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
            % Compute temperature correction factors
            %
            % alpha = obj.TEMPERATURECORRECTIONFACTOR(Tr)
            %
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % alpha : Temperature correction factors for attraction
            % parameters
            omega = obj.AcentricFactor;
            m = 0.3796 + 1.485*omega - 0.1644*omega.^2 + 0.01667*omega.^3;
            alpha = (1 + m.*(1 - sqrt(Tr))).^2;
        end
    end
end