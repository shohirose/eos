classdef VanDerWaalsEos < eos.CubicEosBase
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
            coeffs = [1, -B - 1, A, -A*B];
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
            lnPhi = z - 1 - log(z - B) - A./z;
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
            phi = exp(eos.VanDerWaalsEos.lnFugacityCoeff(z,A,B));
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
            obj@eos.CubicEosBase(0.421875,0.125,Pc,Tc,Mw);
        end
        function obj = setParams(obj,Pc,Tc,Mw)
            % Set parameters
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw);
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
            P = R*T./(V - b) - a./V^2;
        end
        function [z,A,B] = zFactors(obj,P,T)
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
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            Pr = obj.reducedPressure(P);
            Tr = obj.reducedTemperature(T);
            A = obj.reducedAttractionParam(Pr,Tr,1);
            B = obj.reducedRepulsionParam(Pr,Tr);
            x = roots(obj.zFactorCubicEq(A,B));
            z = x(imag(x) == 0);
        end
    end
end