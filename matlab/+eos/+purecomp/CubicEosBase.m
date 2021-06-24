classdef CubicEosBase
    % Base class for two-parameter cubic equation of state.
    
    properties (SetAccess = private)
        CriticalPressure    % Critical pressure [Pa]
        CriticalTemperature % Critical temperature [K]
        MolecularWeight     % Molecular weight [g/mol]
        OmegaA              % Coefficient for attraction parameter
        OmegaB              % Coefficient for repulsion parameter
        AttractionParam     % Attraction parameter
        RepulsionParam      % Repulsion parameter
    end
    methods
        function obj = CubicEosBase(OmegaA,OmegaB,Pc,Tc,Mw)
            % Construct cubic EOS
            %
            % Parameters
            % ----------
            % OmegaA : Coefficient for attraction parameter
            % OmegaB : Coefficient for repulsion parameter
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            arguments
                OmegaA (1,1) {mustBeNumeric}
                OmegaB (1,1) {mustBeNumeric}
                Pc (1,1) {mustBeNumeric}
                Tc (1,1) {mustBeNumeric}
                Mw (1,1) {mustBeNumeric}
            end
            obj.OmegaA = OmegaA;
            obj.OmegaB = OmegaB;
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = OmegaA*(R*Tc)^2/Pc;
            obj.RepulsionParam = OmegaB*R*Tc/Pc;
        end
        function obj = setParams(obj,Pc,Tc,Mw)
            % Set parameters
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            arguments
                obj
                Pc (1,1) {mustBeNumeric}
                Tc (1,1) {mustBeNumeric}
                Mw (1,1) {mustBeNumeric}
            end
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = obj.OmegaA*(R*Tc)^2/Pc;
            obj.RepulsionParam = obj.OmegaB*R*Tc/Pc;
        end
        function Pr = reducedPressure(obj,P)
            % Compute reduced pressure
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            %
            % Returns
            % -------
            % Pr : Reduced pressure
            Pr = P./obj.CriticalPressure;
        end
        function Tr = reducedTemperature(obj,T)
            % Compute reduced temperature
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            %
            % Returns
            % -------
            % Tr : Reduced temperature
            Tr = T./obj.CriticalTemperature;
        end
        function rho = massDensity(obj,P,T,z)
            % Calculate mass density
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            % z : Z-factor
            %
            % Returns
            % -------
            % rho : Mass density [kg/m3]
            R = eos.ThermodynamicConstants.Gas;
            Mw = obj.MolecularWeight;
            rho = P.*Mw*1e-3./(z*R.*T);
        end
        function A = reducedAttractionParam(obj,Pr,Tr,alpha)
            % Compute reduced attraction parameter
            %
            % Parameters
            % ----------
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            % alpha : Temperature correction factor
            %
            % Returns
            % -------
            % A : Reduced attraction parameter
            A = obj.OmegaA*alpha.*Pr./Tr.^2;
        end
        function B = reducedRepulsionParam(obj,Pr,Tr)
            % Compute reduced repulsion parameter
            %
            % Parameters
            % ----------
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % B : Reduced repulsion parameter
            B = obj.OmegaB*Pr./Tr;
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
            arguments
                obj {mustBeA(obj,'eos.purecomp.CubicEosBase')}
                T (:,:) {mustBeNumeric}
                V (:,:) {mustBeNumeric}
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            a = alpha*obj.AttractionParam;
            b = obj.RepulsionParam;
            P = obj.pressureImpl(T,V,a,b);
        end
        function lnPhi = lnFugacityCoeff(obj,z,s)
            % Compute fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factors
            % s : struct containing parameters
            %
            % Returns
            % -------
            % lnPhi : Fugacity coefficients
            arguments
                obj {mustBeA(obj,'eos.purecomp.CubicEosBase')}
                z (:,1) {mustBeNumeric}
                s struct
            end
            lnPhi = exp(obj.lnFugacityCoeffImpl(z,s.A,s.B));
        end
        function phi = fugacityCoeff(obj,z,s)
            % Compute fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factors
            % s : struct containing parameters
            %
            % Returns
            % -------
            % phi : Fugacity coefficients
            arguments
                obj {mustBeA(obj,'eos.purecomp.CubicEosBase')}
                z (:,1) {mustBeNumeric}
                s struct
            end
            phi = exp(obj.lnFugacityCoeffImpl(z,s.A,s.B));
        end
        function [z,s] = zFactors(obj,P,T)
            % Compute Z-factors
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            %
            % Returns
            % -------
            % z : Z-factors
            % s : struct containing parameters
            arguments
                obj {mustBeA(obj,'eos.purecomp.CubicEosBase')}
                P (1,1) {mustBeNumeric}
                T (1,1) {mustBeNumeric}
            end
            Pr = obj.reducedPressure(P);
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            A = obj.reducedAttractionParam(Pr,Tr,alpha);
            B = obj.reducedRepulsionParam(Pr,Tr);
            x = roots(obj.zFactorCubicEq(A,B));
            z = x(imag(x) == 0);
            if nargout > 1
                s.P = P;
                s.T = T;
                s.A = A;
                s.B = B;
            end
        end
        function P = tripleRootPressureRange(obj,T)
            % Compute pressure range with triple roots of Z-factors at a
            % given temperature
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            %
            % Returns
            % -------
            % P : Pressures [Pa]
            arguments
                obj {mustBeA(obj,'eos.purecomp.CubicEosBase')}
                T (1,1) {mustBeNumeric}
            end
            if T >= obj.CriticalTemperature
                error("Error. \nT %f must be less than Tc %f.", ...
                    T, obj.CriticalTemperature);
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            a = alpha*obj.AttractionParam;
            b = obj.RepulsionParam;
            x = roots(obj.dPdVPolyEq(T,a,b));
            V = x(imag(x) == 0);
            V = V(V > b);
            V = sort(V);
            P = obj.pressure(T,V);
        end
    end
end