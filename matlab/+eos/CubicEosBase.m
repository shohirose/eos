classdef CubicEosBase
    % CubicEosBase Base class for two-parameter cubic equation of state.
    %
    %  This class can be used as a basis of two-parameter cubic equation of
    %  state.
    
    properties (SetAccess = private)
        CriticalPressure    % Critical pressure [Pa]
        CriticalTemperature % Critical temperature [K]
        MolecularWeight     % Molecular weight [g/mol]
        OmegaA              % Coefficient for attraction parameter
        OmegaB              % Coefficient for repulsion parameter
        AttractionParam     % Attraction parameter
        RepulsionParam      % Repulsion parameter
    end
    methods (Static)
        function rho = molarDensity(P,T,z)
            % Calculates molar density
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            % z : Z-factor
            %
            % Returns
            % -------
            % rho : Molar density [mol/m3]
            R = eos.ThermodynamicConstants.Gas;
            rho = P/(z*R*T);
        end
        function vm = molarVolume(P,T,z)
            % Calculates molar volume
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            % z : Z-factor
            %
            % Returns
            % -------
            % vm : Molar volume [m3/mol]
            R = eos.ThermodynamicConstants.Gas;
            vm = z*R*T/P;
        end
    end
    methods
        function obj = CubicEosBase(OmegaA,OmegaB,Pc,Tc,Mw)
            % Constructs cubic EOS
            %
            % Parameters
            % ----------
            % OmegaA : Coefficient for attraction parameter
            % OmegaB : Coefficient for repulsion parameter
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
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
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = obj.OmegaA*(R*Tc)^2/Pc;
            obj.RepulsionParam = obj.OmegaB*R*Tc/Pc;
        end
        function Pr = reducedPressure(obj,P)
            % Computes reduced pressure
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            %
            % Returns
            % -------
            % Pr : Reduced pressure
            Pr = P/obj.CriticalPressure;
        end
        function Tr = reducedTemperature(obj,T)
            % Computes reduced temperature
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            %
            % Returns
            % -------
            % Tr : Reduced temperature
            Tr = T/obj.CriticalTemperature;
        end
        function rho = massDensity(obj,P,T,z)
            % Calculates mass density
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
            rho = P*obj.MolecularWeight*1e-3/(z*R*T);
        end
        function A = reducedAttractionParam(obj,Pr,Tr,alpha)
            % Computes reduced attraction parameter
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
            A = obj.OmegaA*alpha*Pr/Tr^2;
        end
        function B = reducedRepulsionParam(obj,Pr,Tr)
            % Computes reduced repulsion parameter
            %
            % Parameters
            % ----------
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % B : Reduced repulsion parameter
            B = obj.OmegaB*Pr/Tr;
        end
    end
end