classdef SoaveRedlichKwongEos < eos.CubicEosBase
    % SoaveRedlichKwongEos Soave-Redlich-Kwong equation of state
    %
    %  This class provides methods to calculate thermodynamic properties
    %  based on Soave-Redlich-Kwong equation of state.
    
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
            coeffs = [1, -1, A - B - B^2, -A*B];
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
            lnPhi = z - 1 - log(z - B) - A/B*log(B./z + 1);
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
            phi = exp(eos.SoaveRedlichKwongEos.lnFugacityCoeff(z,A,B));
        end
    end
    methods
        function obj = SoaveRedlichKwongEos(Pc,Tc,omega,Mw)
            % Constructs SRK EOS
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            obj@eos.CubicEosBase(0.42748,0.08664,Pc,Tc,Mw)
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw)
            % Set parameters
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
            % obj : SoaveRedlichKwongEos
            obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw);
            obj.AcentricFactor = omega;
        end
        function alpha = temperatureCorrectionFactor(obj,Tr)
            % Computes temperature correction factor for attraction parameter
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
            alpha = (1 + m*(1 - sqrt(Tr)))^2;
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
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            a = obj.AttractionParam;
            b = obj.RepulsionParam;
            R = eos.ThermodynamicConstants.Gas;
            P = R*T./(V - b) - alpha*a./(V.*(V + b));
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
            alpha = obj.temperatureCorrectionFactor(Tr);
            A = obj.reducedAttractionParam(Pr,Tr,alpha);
            B = obj.reducedRepulsionParam(Pr,Tr);
            x = roots(eos.SoaveRedlichKwongEos.zFactorCubicEq(A,B));
            z = x(imag(x) == 0);
        end
    end
end