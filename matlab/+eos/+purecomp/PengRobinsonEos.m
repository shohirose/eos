classdef PengRobinsonEos < eos.purecomp.CubicEosBase
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
            % Computes coefficients of Z-factor cubic equation
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
        function coeffs = dPdTPolyEq(T,a,b)
            % Compute coefficients of the polynomial of dPdT = 0.
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
            coeffs = [R*T, 4*b*R*T - 2*a, 2*(b^2*R*T + a*b), ...
                2*b^2*(a - 2*b*R*T), b^3*(b*R*T - 2*a)];
        end
        function lnPhi = lnFugacityCoeff(z,s)
            % Compute the natural log of fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factor
            % s : struct containing parameters
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            arguments
                z (:,1) {mustBeNumeric}
                s struct
            end
            Sqrt2 = eos.purecomp.PengRobinsonEos.Sqrt2;
            Delta1 = eos.purecomp.PengRobinsonEos.Delta1;
            Delta2 = eos.purecomp.PengRobinsonEos.Delta2;
            A = s.A;
            B = s.B;
            lnPhi = z - 1 - log(z - B) ...
                - A./(2*Sqrt2*B).*log((z + Delta1*B)./(z + Delta2*B));
        end
    end
    methods
        function obj = PengRobinsonEos(Pc,Tc,omega,Mw)
            % Construct PR EoS
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
            % obj : PengRobinsonEos
            arguments
                Pc (1,1) {mustBeNumeric}
                Tc (1,1) {mustBeNumeric}
                omega (1,1) {mustBeNumeric}
                Mw (1,1) {mustBeNumeric}
            end
            obj@eos.purecomp.CubicEosBase(0.45724,0.07780,Pc,Tc,Mw);
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw)
            % Set parameters.
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
            % obj : PengRobinsonEos
            arguments
                obj {mustBeA(obj,'eos.purecomp.PengRobinsonEos')}
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
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % alpha : Temperature correction factor
            arguments
                obj {mustBeA(obj,'eos.purecomp.PengRobinsonEos')}
                Tr (:,:) {mustBeNumeric}
            end
            omega = obj.AcentricFactor;
            m = 0.3796 + 1.485*omega - 0.1644*omega^2 + 0.01667*omega^3;
            alpha = (1 + m*(1 - sqrt(Tr))).^2;
        end
        function P = pressure(obj,T,V)
            % Computes pressure
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % V : Volume [m3]
            %
            % Returns
            % -------
            % P : Pressure [Pa]
            arguments
                obj {mustBeA(obj,'eos.purecomp.PengRobinsonEos')}
                T (:,:) {mustBeNumeric}
                V (:,:) {mustBeNumeric}
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            a = obj.AttractionParam;
            b = obj.RepulsionParam;
            R = eos.ThermodynamicConstants.Gas;
            P = R*T./(V - b) - alpha*a./((V - b).*(V + b) + 2*b*V);
        end
    end
end