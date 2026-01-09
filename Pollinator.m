classdef Pollinator < Species
    %POLLINATOR represents a pollinators
    %equivalent to species superclass but with added properties to store
    %(minimal) extincion conditions
    properties
        extinctionConditions = {};
        minimalExtinctionConditions = {};
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
            speciesArray = createArray(size, class(Pollinator));

            for i = 1:numel(speciesArray)
                speciesArray(i) = Pollinator(widthDistribution, constrained);
            end
        end
    end
end

