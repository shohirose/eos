classdef VanDerWaalsEos < eos.purecomp.CubicEosBase
    % VanDerWaalsEos Van der Waals equation of state
    %
    %  This class provides methods to calculate thermodynamic properties
    %  based on Van der Waals equation of state.
    
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
            coeffs = [1, -B - 1, A, -A*B];
        end
        function coeffs = dPdTPolyEq(T,a,b)
            % Compute the coefficients of a polynomial equation of dPdT = 0
            %
            % Parameters
            % ----------
            % T : Temperature
            % a : Attraction parameter
            % b : Repulsion parameter
            %
            % Returns
            % -------
            % coeffs : Coefficients of a polynomial equation
            arguments
                T (1,1) {mustBeNumeric}
                a (1,1) {mustBeNumeric}
                b (1,1) {mustBeNumeric}
            end
            R = eos.ThermodynamicConstants.Gas;
            coeffs = [R*T, -2*a, 4*a*b, -2*a*b^2];
        end
        function lnPhi = lnFugacityCoeff(z,s)
            % Parameters
            % ----------
            % z : Z-factors
            % s : State
            arguments
                z (:,1) {mustBeNumeric}
                s struct
            end
            A = s.A;
            B = s.B;
            lnPhi = z - 1 - log(z - B) - A./z;
        end
    end
    methods
        function obj = VanDerWaalsEos(Pc,Tc,Mw)
            % Constructs VDW EOS
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            % K  : Binary interaction parameters (optional)
            arguments
               Pc (1,1) {mustBeNumeric}
               Tc (1,1) {mustBeNumeric}
               Mw (1,1) {mustBeNumeric}
            end
            obj@eos.purecomp.CubicEosBase(0.421875,0.125,Pc,Tc,Mw);
        end
        function obj = setParams(obj,Pc,Tc,Mw)
            % Set parameters
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            obj = setParams@eos.purecomp.CubicEosBase(obj,Pc,Tc,Mw);
        end
        function alpha = temperatureCorrectionFactor(~,~)
            % Computes temperature correction factor for attraction parameter
            %
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % alpha : Temperature correction factor
            alpha = 1;
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
            R = eos.ThermodynamicConstants.Gas;
            a = obj.AttractionParam;
            b = obj.RepulsionParam;
            P = R*T./(V - b) - a./V.^2;
        end
        %{
        function [z,s] = zFactors(obj,P,T)
            % Computes Z-factors
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            %
            % Returns
            % -------
            % z : Z-factors
            % s : State
            Pr = obj.reducedPressure(P);
            Tr = obj.reducedTemperature(T);
            A = obj.reducedAttractionParam(Pr,Tr,1);
            B = obj.reducedRepulsionParam(Pr,Tr);
            x = roots(obj.zFactorCubicEq(A,B));
            z = x(imag(x) == 0);
            s = struct('A',A,'B',B);
        end
        %}
        %{
        function P = tripleRootPressureRange(obj,T)
            % Computes pressure range with triple roots of Z-factors at a
            % given temperature
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            %
            % Returns
            % -------
            % P : Pressures [Pa]
            if T >= obj.CriticalTemperature
                error("Error. \nT %f must be less than Tc %f.", ...
                    T, obj.CriticalTemperature);
            end
            a = obj.AttractionParam;
            b = obj.RepulsionParam;
            R = eos.ThermodynamicConstants.Gas;
            x = roots([R*T, -2*a, 4*a*b, -2*a*b^2]);
            V = x(imag(x) == 0);
            V = V(V > b);
            V = sort(V);
            P = obj.pressure(T,V);
        end
        %}
    end
end