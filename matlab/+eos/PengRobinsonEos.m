classdef PengRobinsonEos < eos.CubicEosBase
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
        function obj = PengRobinsonEos(Pc,Tc,omega,Mw)
            % Constructs PR EOS
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            obj@eos.CubicEosBase(0.45724,0.07780,Pc,Tc,Mw);
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
            m = 0.3796 + 1.485*omega - 0.1644*omega^2 + 0.01667*omega^3;
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
            P = R*T./(V - b) - alpha*a./((V - b).*(V + b) + 2*b*V);
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
            x = roots(eos.PengRobinsonEos.zFactorCubicEq(A,B));
            z = x(imag(x) == 0);
        end
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
                error("Error. \nTemperature %f must be greater than critical temperature %f.", T, obj.CriticalTemperature);
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            a = alpha*obj.AttractionParam;
            b = obj.RepulsionParam;
            R = eos.ThermodynamicConstants.Gas;
            x = roots([R*T, 2*(2*b*R*T - a), 2*((2*b - 1)*b*R*T + a*b), 2*b^2*(a - 2*R*T), b^2*(R*T - a*b)]);
            V = x(imag(x) == 0);
            V = V(V > b);
            V = sort(V);
            P = obj.pressure(T,V);
        end
    end
end