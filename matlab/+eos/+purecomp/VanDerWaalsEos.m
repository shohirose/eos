classdef VanDerWaalsEos < eos.purecomp.CubicEosBase
    % Van der Waals equation of state.
    
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
            coeffs = [1, -B - 1, A, -A*B];
        end
        function coeffs = dPdTPolyEq(T,a,b)
            % Compute the coefficients of a polynomial equation of dPdT = 0
            %
            % coeffs = DPDTPOLYEQ(T,a,b)
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
            % Compute the natural log of fugacity coeffcients
            %
            % lnPhi = LNGUGACITYCOEFF(z,s)
            %
            % Parameters
            % ----------
            % z : Z-factors
            % s : struct containing parameters
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
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
            % Construct VDW EoS
            %
            % obj = VANDERWAALSEOS(Pc,Tc,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            % K  : Binary interaction parameters (optional)
            %
            % Returns
            % -------
            % obj : VANDERWAALSEOS
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
            % obj = obj.SETPARAMS(Pc,Tc,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            %
            % Returns
            % -------
            % obj : VANDERWAALSEOS
            obj = setParams@eos.purecomp.CubicEosBase(obj,Pc,Tc,Mw);
        end
        function alpha = temperatureCorrectionFactor(~,~)
            % Compute temperature correction factor
            %
            %   This function just returns 1 because VDW EoS does not
            %   consider temperature dependence of attraction parameter.
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
            alpha = 1;
        end
        function P = pressure(obj,T,V)
            % Compute pressure
            %
            % P = obj.PRESSURE(T,V)
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
    end
end