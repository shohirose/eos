classdef CubicEosBase
    % CubicEosBase Base class for two-parameter cubic equation of state.
    %
    %  This class can be used as a basis of two-parameter cubic equation of
    %  state.
    
    properties (SetAccess = private)
        CriticalPressure    % Critical pressure
        CriticalTemperature % Critical temperature
        OmegaA              % Coefficient for attraction parameter
        OmegaB              % Coefficient for repulsion parameter
        AttractionParam     % Attraction parameter
        RepulsionParam      % Repulsion parameter
    end
    methods
        function obj = CubicEosBase(OmegaA,OmegaB,Pc,Tc)
            % Constructs cubic EOS
            %
            % Parameters
            % ----------
            % OmegaA : Coefficient for attraction parameter
            % OmegaB : Coefficient for repulsion parameter
            % Pc : Critical pressure
            % Tc : Critical temperature
            obj.OmegaA = OmegaA;
            obj.OmegaB = OmegaB;
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = OmegaA*(R*Tc)^2/Pc;
            obj.RepulsionParam = OmegaB*R*Tc/Pc;
        end
        function obj = setCriticalProperties(obj,Pc,Tc)
            % Set critical pressure and temperature
            %
            % Parameters
            % ----------
            % Pc : Critical pressure
            % Tc : Critical temperature
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = obj.OmegaA*(R*Tc)^2/Pc;
            obj.RepulsionParam = obj.OmegaB*R*Tc/Pc;
        end
        function Pr = reducedPressure(obj,P)
            % Computes reduced pressure
            %
            % Parameters
            % ----------
            % P : Pressure
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
            % T : Temperature
            %
            % Returns
            % -------
            % Tr : Reduced temperature
            Tr = T/obj.CriticalTemperature;
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