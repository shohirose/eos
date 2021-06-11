classdef CubicEosBase
    properties (SetAccess = private)
        CriticalPressure
        CriticalTemperature
        OmegaA          % Coefficient for attraction parameter
        OmegaB          % Coefficient for repulsion parameter
        AttractionParam
        RepulsionParam
    end
    methods
        function obj = CubicEosBase(OmegaA,OmegaB,Pc,Tc)
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
            Pr = P/obj.CriticalPressure;
        end
        function Tr = reducedTemperature(obj,T)
            % Computes reduced temperature
            %
            % Parameters
            % ----------
            % T : Temperature
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
            A = obj.OmegaA*alpha*Pr/Tr^2;
        end
        function B = reducedRepulsionParam(obj,Pr,Tr)
            % Computes reduced repulsion parameter
            %
            % Parameters
            % ----------
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            B = obj.OmegaB*Pr/Tr;
        end
    end
end