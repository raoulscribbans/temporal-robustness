classdef Interaction
    %INTERACTION represents an interaction between a plant and pollinator    
    properties 
        startTime
        endTime
        plantIndex
        pollinatorIndex
    end
    
    methods
        function obj = Interaction(startTime, endTime, plantIndex, pollinatorIndex)
            %INTERACTION Construct an instance of this class
            arguments
                startTime (1, 1) double = 0;
                endTime (1, 1) double {mustBeGreaterThanOrEqual(endTime, startTime)} = 0
                plantIndex (1, 1) uint16 = 0
                pollinatorIndex (1, 1) uint16 = 0 
            end

            obj.startTime = startTime;
            obj.endTime = endTime;
            obj.plantIndex = plantIndex;
            obj.pollinatorIndex = pollinatorIndex;
        end
    end
end

