import ThermodynamicConstants.Gas

classdef CubicEosBase
    properties
        CriticalPressure
        CriticalTemperature
    end
    properties (Access = protected)
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
            obj.OmegaA = OmegaA
            obj.OmegaB = OmegaB
            obj.CriticalPressure = Pc
            obj.CriticalTemperature = Tc
            obj.AttractionParam = OmegaA*(Gas*Tc)^2/Pc
            obj.RepulsionParam = OmegaB*Gas*Tc/Pc
        end
        function Pr = reducedPressure(obj,P)
            % P : Pressure
            Pr = P/obj.CriticalPressure
        end
        function Tr = reducedTemperature(obj,T)
            % T : Temperature
            Tr = T/obj.CriticalTemperature
        end
        function A = reducedAttractionParam(obj,Pr,Tr,alpha)
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            % alpha : Temperature correction factor
            A = obj.OmegaA*alpha*Pr/Tr^2
        end
        function B = reducedRepulsionParam(obj,Pr,Tr)
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            B = obj.OmegaB*Pr/Tr
        end
    end
end