classdef Species
    %SPECIES Class representing a general species
    %   Species represents either a plant or pollinator, however
    %   pollinators should be represented with the Pollinator subclass

    properties
        time (1,1) double;
        intervalWidth (1,1) double {mustBeNonnegative};
        degree (1, 1) uint16 = 0;
        index (1, 1) uint16 = 0;
        order (1, 1) string = "";
        family (1, 1) string = "";
        genus (1, 1) string = "";
        name (1, 1) string = "";
    end
    
    properties (Dependent)
        startTime;
        endTime;
    end

    methods
        function obj = Species(widthDistribution, constrained, order, family, genus, name)
            %SPECIES construct an instance of this class
            %   randomly assigns an interval width based on the provided
            %   witdth distribution. If constrained is true then the entire
            %   interval will be placed inside [0, 1], but if not the
            %   interval midpoint will be instead.
            arguments
                widthDistribution (1, 1) prob.ProbabilityDistribution = makedist("Normal", mu = 0.1, sigma = 0);
                constrained (1, 1) logical = true
                order (1, 1) string = "";
                family (1, 1) string = "";
                genus (1, 1) string = "";
                name (1, 1) string = "";
            end
            
            obj.intervalWidth = widthDistribution.random;

            if constrained
                obj.time = (obj.intervalWidth / 2) + (1 - obj.intervalWidth) * rand();
            else
                obj.time = rand();
            end

            obj.order = order;
            obj.family = family;
            obj.genus = genus;
            obj.name = name;
        end

        function startTime = get.startTime(obj)
            startTime = obj.time - obj.intervalWidth/2;
        end
        function endTime = get.endTime(obj)
            endTime = obj.time + obj.intervalWidth/2;
        end
    end

    methods (Static)
        function speciesArray = createSpeciesArray(widthDistribution, constrained, size)
            %CREATESPECIESARRAY generate an array (of specified size) of species
            %   species are assigned same width distribution and
            %   constrained value
            arguments
                widthDistribution (1, 1) prob.ProbabilityDistribution = makedist("Normal", mu = 0.1, sigma = 0);
                constrained (1, 1) logical = false;
                size int16 {mustBeVector} = 1;
            end
            speciesArray = createArray(size, class(Species));

            for i = 1:numel(speciesArray)
                speciesArray(i) = Species(widthDistribution, constrained);
            end
        end
    end
end

