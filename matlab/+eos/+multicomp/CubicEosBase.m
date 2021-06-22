classdef CubicEosBase
    % Base class for two-parameter cubic equations of state for
    % multi-component systems.
    
    properties (SetAccess = private)
        CriticalPressure            % Critical pressure [Pa]
        CriticalTemperature         % Critical temperature [K]
        MolecularWeight             % Molecular weight [g/mol]
        OmegaA                      % Coefficient for attraction parameter
        OmegaB                      % Coefficient for repulsion parameter
        AttractionParam             % Attraction parameter
        RepulsionParam              % Repulsion parameter
        BinaryInteractionParams     % Binary interaction parameters
    end
    
    methods
        function obj = CubicEosBase(OmegaA,OmegaB,Pc,Tc,Mw,K)
            % Construct cubic EOS.
            %
            % obj = CUBICEOSBASE(OmegaA,OmegaB,Pc,Tc,Mw,K)
            %
            % Parameters
            % ----------
            % OmegaA : Coefficient for attraction parameter
            % OmegaB : Coefficient for repulsion parameter
            % Pc     : Critical pressure [Pa]
            % Tc     : Critical temperature [K]
            % Mw     : Molecular weight [g/mol]
            % K      : Binary interaction parameters
            arguments
                OmegaA (1,1) {mustBeNumeric}
                OmegaB (1,1) {mustBeNumeric}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric}
            end
            obj.OmegaA = OmegaA;
            obj.OmegaB = OmegaB;
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = OmegaA*R^2*Tc.^2./Pc;
            obj.RepulsionParam = OmegaB*R*Tc./Pc;
            obj.BinaryInteractionParams = K;
        end
        function obj = setParams(obj,Pc,Tc,Mw,K)
            % Set parameters.
            %
            % obj = obj.SETPARAMS(Pc,Tc,Mw,K)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            % K  : Binary interaction parameters
            %
            % Returns
            % -------
            % obj : an instance of CubicEosBase
            arguments
                obj {mustBeA(obj,'eos.multicomp.CubicEosBase')}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric}
            end
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = obj.OmegaA*R^2*Tc.^2./Pc;
            obj.RepulsionParam = obj.OmegaB*R*Tc./Pc;
            obj.BinaryInteractionParams = K;
        end
        function Pr = reducedPressure(obj,P)
            % Compute reduced pressure.
            %
            % Pr = obj.REDUCEDPRESSURE(P)
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
            % Compute reduced temperature.
            %
            % Tr = obj.REDUCEDTEMPERATURE(T)
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
        function rho = massDensity(obj,P,T,z,x)
            % Calculate mass density.
            %
            % rho = obj.MASSDENSITY(P,T,z,x)
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            % z : Z-factor
            % x : Composition
            %
            % Returns
            % -------
            % rho : Mass density [kg/m3]
            R = eos.ThermodynamicConstants.Gas;
            Mw = x'*obj.MolecularWeight;
            rho = P.*Mw*1e-3./(z*R.*T);
        end
        function A = reducedAttractionParam(obj,Pr,Tr,alpha)
            % Compute reduced attraction parameter.
            %
            % A = obj.REDUCEDATTRACTIONPARAM(Pr,Tr,alpha)
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
            % Computes reduced repulsion parameter.
            %
            % B = obj.REDUCEDREPULSIONPARAM(Pr,Tr)
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
        function P = pressure(obj,T,V,x)
            % Compute pressure.
            %
            % P = obj.PRESSURE(T,V,x)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % V : Volume [m3]
            % x : Phase composition
            %
            % Returns
            % -------
            % P : Pressure [Pa]
            arguments
                obj {mustBeA(obj,'eos.multicomp.CubicEosBase')}
                T (1,1) {mustBeNumeric}
                V (:,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric} 
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            ai = alpha.*obj.AttractionParam;
            bi = obj.RepulsionParam;
            [a,b] = obj.applyMixingRule(x,ai,bi);
            P = obj.pressureImpl(T,V,a,b);
        end
        function [z,s] = zFactors(obj,P,T,x)
            % Computes Z-factors.
            %
            % [z,s] = obj.ZFACTORS(P,T,x)
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            % x : Composition
            %
            % Returns
            % -------
            % z : Z-factors
            % s : struct containing parameters
            Pr = obj.reducedPressure(P);
            Tr = obj.reducedTemperature(T);
            Ai = obj.reducedAttractionParam(Pr,Tr,1);
            Bi = obj.reducedRepulsionParam(Pr,Tr);
            [A,B,Aij] = obj.applyMixingRule(x,Ai,Bi);
            y = roots(obj.zFactorCubicEq(A,B));
            z = y(imag(y) == 0);
            if nargout > 1
                s.P = P;
                s.T = T;
                s.x = x;
                s.A = A;
                s.B = B;
                s.Ai = Ai;
                s.Bi = Bi;
                s.Aij = Aij;
            end
        end
        function phi = fugacityCoeff(obj,z,s)
            % Compute fugacity coefficients.
            %
            % phi = obj.FUGACITYCOEFF(z,s)
            %
            % Parameters
            % ----------
            % z : Z-factors
            % s : struct containg parameters
            %
            % Returns
            % -------
            % phi : Fugacity coefficients
            arguments
                obj {mustBeA(obj,'eos.multicomp.CubicEosBase')}
                z (:,1) {mustBeNumeric}
                s struct
            end
            phi = exp(obj.lnFugacityCoeff(z,s));
        end
        function [a,b,aij] = applyMixingRule(obj,x,ai,bi)
            % Apply mixing rule to attraction and repulsion parameters.
            %
            % [a,b,aij] = obj.APPLYMIXINGRULE(x,ai,bi)
            %
            % Parameters
            % ----------
            % x  : Composition
            % ai : Attraction parameter of each component
            % bi : Repulsion parameter of each component
            %
            % Returns
            % -------
            % a   : Attraction parameter of the mixture
            % b   : Repulsion parameter of the mixture
            % aij : Combined attraction parameter between i and j
            % components
            K = obj.BinaryInteractionParams;
            % Combining rule with correction parameters
            aij = (1 - K).*sqrt(kron(ai,ai'));
            % Quadratic mixing
            a = x'*aij*x;
            % Linear mixing
            b = x'*bi;
        end
        function P = tripleRootPressureRange(obj,T,x)
            % Compute a pressure range with triple roots of Z-factors
            %
            % P = obj.TRIPLEROOTPRESSURERANGE(T,x)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % x : Phase composition
            %
            % Returns
            % -------
            % P : Pressures [Pa]
            arguments
                obj {mustBeA(obj,'eos.multicomp.CubicEosBase')}
                T (1,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric}
            end
            if T >= obj.CriticalTemperature
                error("Error. \nT %f must be less than Tc %f.", ...
                    T, obj.CriticalTemperature);
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            ai = alpha.*obj.AttractionParam;
            bi = obj.RepulsionParam;
            [a,b] = obj.applyMixingRule(x,ai,bi);
            y = roots(obj.dPdTPolyEq(T,a,b));
            V = y(imag(y) == 0);
            V = V(V > b);
            V = sort(V);
            P = obj.pressure(T,V,x);
        end
    end
    methods (Abstract)
        alpha = temperatureCorrectionFactor(obj,Tr)
    end
end