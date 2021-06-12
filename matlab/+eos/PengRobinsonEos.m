classdef PengRobinsonEos < eos.CubicEosBase
    properties (Constant)
        Sqrt2 = sqrt(2)
        Delta1 = 1 + sqrt(2)
        Delta2 = 1 - sqrt(2)
    end
    properties
        AcentricFactor
    end
    methods (Static)
        function coeffs = zFactorCubicEq(A,B)
            % Computes coefficients of Z-factor cubic equation
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            coeffs = [1, B - 1, A - 2*B - 3*B^2, -A*B + B^2 + B^3];
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
            Sqrt2 = eos.PengRobinsonEos.Sqrt2;
            Delta1 = eos.PengRobinsonEos.Delta1;
            Delta2 = eos.PengRobinsonEos.Delta2;
            lnPhi = z - 1 - log(z - B) - A./(2*Sqrt2*B).*log((z + Delta1*B)./(z + Delta2*B));
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
            phi = exp(eos.PengRobinsonEos.lnFugacityCoeff(z,A,B));
        end
    end
    methods
        function obj = PengRobinsonEos(Pc,Tc,omega)
            % Pc : Critical pressure
            % Tc : Critical temperature
            % omega : Acentric factor
            obj@eos.CubicEosBase(0.45724,0.07780,Pc,Tc);
            obj.AcentricFactor = omega;
        end
        function obj = setCriticalProperties(obj,Pc,Tc,omega)
            % Set critical properties
            %
            % Parameters
            % ----------
            % Pc : Critical pressure
            % Tc : Critical temperature
            % omega : Acentric factor
            obj = setCriticalProperties@eos.CubicEosBase(obj,Pc,Tc);
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
            m = 0.3796 + 1.485*omega - 0.1644*omega^2 + 0.01667*omega^3;
            alpha = (1 + m*(1 - sqrt(Tr)))^2;
        end
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
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            a = obj.AttractionParam;
            b = obj.RepulsionParam;
            R = eos.ThermodynamicConstants.Gas;
            P = R*T./(V - b) - alpha*a./((V - b).*(V + b) + 2*b*V);
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
            Pr = obj.reducedPressure(P);
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            A = obj.reducedAttractionParam(Pr,Tr,alpha);
            B = obj.reducedRepulsionParam(Pr,Tr);
            x = roots(eos.PengRobinsonEos.zFactorCubicEq(A,B));
            z = x(imag(x) == 0);
        end
    end
end