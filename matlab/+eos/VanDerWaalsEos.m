import CubicEosBase
import ThermodynamicConstants.Gas

classdef VanDerWaalsEos < CubicEosBase
    methods
        function obj = VanDerWaalsEos(obj,Pc,Tc)
            % Pc : Critical pressure
            % Tc : Critical temperature
            obj@CubicEosBase(0.421875,0.125,Pc,Tc)
        end
    end
    methods (Static)
        function coeffs = zFactorCubicEq(A,B)
            % Computes coefficients of Z-factor cubic equation
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            coeffs = [1, -B - 1, A, -A*B]
        end
        function lnPhi = lnFugacityCoeff(z,A,B)
            % Computes natural log of fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factor
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            lnPhi = z - 1 - log(z - B) - A./z
        end
        function phi = fugacityCoeff(z,A,B)
            % Computes fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factor
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            %
            % Returns
            % -------
            % phi : Fugacity coefficients
            phi = exp(lnFugacityCoeff(z,A,B))
        end
    end
    methods
        function P = pressure(obj,T,V)
            % Computes pressure
            %
            % Parameters
            % ----------
            % T : Temperature
            % V : Volume
            %
            % Returns
            % -------
            % P : Pressure
            a = obj.attractionParam
            b = obj.repulsionParam
            P = Gas*T./(V - b) - a./V^2
        end
        function [z,A,B] = zFactors(obj,P,T)
            % Computes Z-factors
            %
            % Parameters
            % ----------
            % P : Pressure
            % T : Temperature
            %
            % Returns
            % -------
            % z : Z-factors
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            Pr = obj.reducedPressure(P)
            Tr = obj.reducedTemperature(T)
            A = obj.reducedAttractionParam(Pr,Tr,1)
            B = obj.reducedRepulsionParam(Pr,Tr)
            x = roots(zFactorCubicEq(A,B))
            z = x(imag(x) == 0)
        end
    end
end