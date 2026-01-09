classdef TemporalNetwork
    %TEMPORALNETWORK contains a set of plants, pollinators, and their
    %interactions, with methods to evaluate network metrics including
    %robustness
    
    properties
        plants Species 
        pollinators Pollinator
        interactions Interaction
        robustness double
        interactionProbability double
        minimumReproductionTime double
        overlapFunction function_handle
    end
    
    properties (Dependent)
        numberOfPlants
        numberOfPollinators
    end

    methods
        function obj = TemporalNetwork(plants, pollinators, interactionProbability, minimumReproductionTime, overlapFunction)
            %TEMPORALNETWORK Construct an instance of this class
            %overlap function determines how interaction probability
            %depends on overlap, and will overrule interaction probability,
            %EXCEPT in the case where it is not specified
            arguments
                plants Species {mustBeVector(plants)}
                pollinators Pollinator {mustBeVector(pollinators)}
                interactionProbability double {mustBeInRange(interactionProbability, 0, 1)} = 1
                minimumReproductionTime double {mustBeInRange(minimumReproductionTime, 0, 1)} = 1
                overlapFunction function_handle = @(x)(x>0)*interactionProbability %default function applies equal p_int for all overlaps
            end

            obj.plants = plants;
            plantIndices = num2cell(1:obj.numberOfPlants);
            [obj.plants.index] = plantIndices{:};

            obj.pollinators = pollinators;
            pollinatorIndices = num2cell(1:obj.numberOfPollinators);
            [obj.pollinators.index] = pollinatorIndices{:};

            obj.interactionProbability = interactionProbability;
            obj.minimumReproductionTime = minimumReproductionTime;
            obj.overlapFunction = overlapFunction;

            obj = obj.setInteractionsAndDegrees();
            obj = obj.setExtinctionConditions();
        end
        
        function numberOfPlants = get.numberOfPlants(obj)
            numberOfPlants = length(obj.plants);
        end
        
        function numberOfPollinators = get.numberOfPollinators(obj)
            numberOfPollinators = length(obj.pollinators);
        end
        
        function adjacencyMatrix = getWeightedAdjacencyMatrix(obj)
            adjacencyMatrix = zeros(obj.numberOfPlants, obj.numberOfPollinators);

            for i = 1:length(obj.interactions)
                adjacencyMatrix(obj.interactions(i).plantIndex, obj.interactions(i).pollinatorIndex) = obj.interactions(i).endTime - obj.interactions(i).startTime;
            end
        end
            
        function adjacencyMatrix = getBinaryAdjacencyMatrix(obj)
            adjacencyMatrix = logical(getWeightedAdjacencyMatrix(obj));
        end
    
        function robustnessData = simulateExtinctions(obj, numberOfIterations)
            
            robustnessData = struct( ...
                'numberOfIterations', numberOfIterations);
            
            meanPollinatorsRemaining = zeros(obj.numberOfPlants + 1, 1);
            robustnesses = zeros(numberOfIterations, 1);

            for it = 1:numberOfIterations

                extinctionSequence = randperm(obj.numberOfPlants);
                pollinatorsRemaining = zeros(obj.numberOfPlants + 1, 1);

                isPollinatorAlive = true(obj.numberOfPollinators, 1);

                for i = 1:obj.numberOfPlants
                    for j = 1:obj.numberOfPollinators
                        if isPollinatorAlive(j) && obj.isPollinatorExtinct( ...
                                extinctionSequence(1:i-1), ...
                                obj.pollinators(j))
                            isPollinatorAlive(j) = false;
                        end
                    end
                    pollinatorsRemaining(i) = sum(isPollinatorAlive);
                end

                meanPollinatorsRemaining = meanPollinatorsRemaining + pollinatorsRemaining / numberOfIterations;
                robustnesses(it) = sum(pollinatorsRemaining)/((obj.numberOfPlants + 1) * obj.numberOfPollinators);
            end
            
            robustnessData.meanPollinatorsRemaining = meanPollinatorsRemaining;
            robustnessData.robustnesses = robustnesses;
        end
    
        function robustness = calculateRobustness(obj)
            robustness = mean(obj.calculateIndividualRobustness());
        end
        
        function individualRobustness = calculateIndividualRobustness(obj)
        
            individualRobustness = ones(obj.numberOfPollinators, 1);

            for i=1:obj.numberOfPollinators

                minimalExtinctionConditions = obj.pollinators(i).minimalExtinctionConditions;

                if isempty(minimalExtinctionConditions) || isempty(minimalExtinctionConditions{1})
                    individualRobustness(i) = 0;
                else
                    numberOfMinimalConditions = length(minimalExtinctionConditions);
                    for l = 1:numberOfMinimalConditions

                        sign = (-1)^(l+1);
                        combinations = nchoosek(1:numberOfMinimalConditions, l);

                        lengths = arrayfun(@(k) length(unique([minimalExtinctionConditions{combinations(k, :)}])), 1:size(combinations, 1));

                        individualRobustness(i) = individualRobustness(i) - sign * sum(1 ./ (1 + lengths));
                    end
                end
            end
        end

        function obj = setInteractions(obj, adjacencyMatrix)
            arguments
                obj
                adjacencyMatrix logical
            end

            plantStartTimes = [obj.plants.startTime];
            plantEndTimes = [obj.plants.endTime];
            pollinatorStartTimes = [obj.pollinators.startTime];
            pollinatorEndTimes = [obj.pollinators.endTime];

            newInteractions(obj.numberOfPlants * obj.numberOfPollinators, 1) = Interaction;
            numberOfNewInteractions = 0;
            removedInteractions = false(1 ,length(obj.interactions));

            for i = 1:obj.numberOfPlants
                for j = 1:obj.numberOfPollinators
                    interactionIndex = ([obj.interactions.plantIndex] == i) & ([obj.interactions.pollinatorIndex] == j);
                    interaction = obj.interactions(interactionIndex);

                    if adjacencyMatrix(i, j) && isempty(interaction)
                        latestStart = max(plantStartTimes(i), pollinatorStartTimes(j));
                        earliestFinish = min(plantEndTimes(i), pollinatorEndTimes(j));

                        if latestStart > earliestFinish
                            warning("impossible interaction between pollinator" + j + " and plant " + i);
                        else
                            numberOfNewInteractions = numberOfNewInteractions + 1;
                            newInteractions(numberOfNewInteractions) = Interaction(latestStart, earliestFinish, i, j);
                        end
                    elseif ~adjacencyMatrix(i, j) && ~isempty(interaction)
                        removedInteractions = removedInteractions | interactionIndex;
                    end
                end
            end

            obj.interactions = [obj.interactions(~removedInteractions); newInteractions(1:numberOfNewInteractions)];
            obj = obj.setExtinctionConditions();
        end
        
        function timedDegrees = getTimedDegrees(obj, timeStep, setInactiveNan)
            arguments
                obj
                timeStep (1, 1) double {mustBeLessThan(timeStep, 1), mustBePositive}
                setInactiveNan (1, 1) logical = false
            end

            times = 0:timeStep:1;
            timedDegrees = zeros(length(times), obj.numberOfPollinators);
            adjacencyMatrix = obj.getBinaryAdjacencyMatrix();

            plantStartTimes = [obj.plants.startTime];
            plantEndTimes = [obj.plants.endTime];

            pollinatorStartTimes = [obj.pollinators.startTime];
            pollinatorEndtimes = [obj.pollinators.endTime];

            for t = 1:length(times)
                plantsActive = plantStartTimes < times(t) & plantEndTimes > times(t);
                pollinatorsActive = pollinatorStartTimes < times(t) & pollinatorEndtimes > times(t);             

                timedDegrees(t, :) = (plantsActive * adjacencyMatrix) .* pollinatorsActive;

                if setInactiveNan
                    timedDegrees(t, ~pollinatorsActive) = NaN;
                end
            end
        end

        function timeAveragedDegrees = getTimeAveragedDegrees(obj, timeStep)
            timeAveragedDegrees = mean(obj.getTimedDegrees(timeStep)) ./ [obj.pollinators.intervalWidth];
        end

        function overlapMatrix = getOverlapMatrix(obj)

            plantStartTimes = [obj.plants.startTime];
            plantEndTimes = [obj.plants.endTime];

            pollinatorStartTimes = [obj.pollinators.startTime];
            pollinatorEndtimes = [obj.pollinators.endTime];

            overlapMatrix = zeros(obj.numberOfPlants, obj.numberOfPollinators);

            for i = 1:obj.numberOfPlants
                for j = 1:obj.numberOfPollinators
                    overlaps = plantStartTimes(i) < pollinatorEndtimes(j) && plantEndTimes(i) > pollinatorStartTimes(j);
                    if overlaps
                        sortedTimes = sort([plantStartTimes(i), pollinatorEndtimes(j), plantEndTimes(i), pollinatorStartTimes(j)]);
                        overlapMatrix(i, j) = sortedTimes(3) - sortedTimes(2); %overlap is given by difference of two centremost times
                    end
                end
            end
        end

        function binaryOverlapMatrix = getBinaryOverlapMatrix(obj)
            binaryOverlapMatrix = obj.getOverlapMatrix > 0;
        end
    
        function nullNetworks = getUniformNullModelNetwork(obj, iterations, pollinatorWidth, breedingFraction, constrained, overlapFunction, plantWidthDistribution)
            arguments
                obj 
                iterations (1, 1) uint16 = 1  
                pollinatorWidth (1, 1) double = mean([obj.pollinators.intervalWidth])
                breedingFraction (1, 1) double = 0
                constrained (1, 1) logical = true;
                overlapFunction (1, 1) function_handle = obj.overlapFunction
                plantWidthDistribution (1, 1) prob.ProbabilityDistribution = makedist("Normal", mu = mean([obj.plants.intervalWidth]), sigma = 0)
            end

            pollinatorWidthDistribution = makedist("Normal", mu = pollinatorWidth, sigma = 0);
            
            estimatedInteractionProbability = length(obj.interactions) / sum(obj.getBinaryOverlapMatrix, 'all');

            for i = 1:iterations
                for j = length(obj.plants):-1:1
                    nullPlants(j) = Species(plantWidthDistribution, constrained);
                end
                for j = length(obj.pollinators):-1:1
                    nullPollinators(j) = Pollinator(pollinatorWidthDistribution, constrained);
                end

                if i == 1
                    nullNetworks = repmat(TemporalNetwork(nullPlants, nullPollinators, estimatedInteractionProbability, breedingFraction, overlapFunction), iterations, 1);
                else
                    nullNetworks(i) = TemporalNetwork(nullPlants, nullPollinators, estimatedInteractionProbability, breedingFraction, overlapFunction);
                end
            end
        end

        function plot = timePlot(obj, barWidth, plantColour, pollinatorColour)
            plot = plotHelper.timePlot(obj.plants, obj.pollinators, barWidth, plantColour, pollinatorColour);
        end

        function p = graphPlot(obj, plantColour, pollinatorColour, name)
            arguments
                obj 
                plantColour 
                pollinatorColour 
                name = true 
            end
            [~, orderPlantIndices] = sort([obj.plants.time], 'descend');
            [~, orderPollinatorsIndices] = sort([obj.pollinators.time], 'ascend');
            
            weightedAdjacency = obj.getWeightedAdjacencyMatrix();
            
            bigAdjacencyMatrix = [zeros(obj.numberOfPlants), weightedAdjacency(orderPlantIndices, orderPollinatorsIndices); weightedAdjacency(orderPlantIndices, orderPollinatorsIndices)', zeros(obj.numberOfPollinators)];
            
            if name
                g = graph(bigAdjacencyMatrix, [obj.plants.name; obj.pollinators.name]);
            else
                g = graph(bigAdjacencyMatrix);
            end

            markerScaleFactor = 10;
            edgeScaleFactor = 5;

            p = plot(g, ...,
                "Layout","circle",...
                "NodeColor", [repmat(plantColour, obj.numberOfPlants, 1); repmat(pollinatorColour, obj.numberOfPollinators, 1)], ...
                'LineWidth', g.Edges.Weight * edgeScaleFactor, ... 
                MarkerSize = [obj.plants(orderPlantIndices).intervalWidth, obj.pollinators(orderPollinatorsIndices).intervalWidth] * markerScaleFactor, Marker="o", ...
                Marker=[repmat("o", 1, obj.numberOfPlants), repmat("square", 1, obj.numberOfPollinators)], ...
                NodeFontSize=10);

            if ~name
                labelnode(p,1:(obj.numberOfPlants + obj.numberOfPollinators),'');
            end
        end
    end

    methods(Static)
        function networks = createFromDistributions( ...
                plantWidthDistribution, ...
                pollinatorWidthDistribution, ...
                numberOfPlants, ...
                numberOfPollinators, ...
                constrained, ...
                size, ...
                interactionProbability, ...
                minimumReproductionTime, ...
                overlapFunction)
            arguments
                plantWidthDistribution 
                pollinatorWidthDistribution 
                numberOfPlants 
                numberOfPollinators 
                constrained (1, 1) logical = true
                size {mustBeVector, mustBeInteger} = 1
                interactionProbability (1, 1) double {mustBeInRange(interactionProbability, 0, 1)} = 1
                minimumReproductionTime double {mustBeInRange(minimumReproductionTime, 0, 1)} = 1
                overlapFunction function_handle = @(x)(x>0)*interactionProbability %default function applies equal p_int for all overlaps
            end

            networks = repmat(TemporalNetwork(Species(), Pollinator()), size);

            for i = 1:numel(networks)
                plants = Species.createSpeciesArray(plantWidthDistribution, constrained, [numberOfPlants, 1]);
                pollinators = Pollinator.createSpeciesArray(pollinatorWidthDistribution, constrained, [numberOfPollinators, 1]);
                networks(i) = TemporalNetwork(plants, pollinators, interactionProbability, minimumReproductionTime, overlapFunction);
            end
        end
    end
    methods (Access = private)
        function obj = setInteractionsAndDegrees(obj)

            interactions(obj.numberOfPlants * obj.numberOfPollinators, 1) = Interaction;
            numberOfInteractions = 0;
            
            plantStartTimes = [obj.plants.startTime];
            plantEndTimes = [obj.plants.endTime];
            pollinatorStartTimes = [obj.pollinators.startTime];
            pollinatorEndTimes = [obj.pollinators.endTime];
            for i = 1:obj.numberOfPlants   
                for j = 1:obj.numberOfPollinators
                    latestStart = max(plantStartTimes(i), pollinatorStartTimes(j));
                    earliestFinish = min(plantEndTimes(i), pollinatorEndTimes(j));
            
                    if(latestStart < earliestFinish && rand <= obj.overlapFunction(earliestFinish - latestStart))
                        if obj.overlapFunction(earliestFinish - latestStart) > 1
                            warning("overlap prob = " + obj.overlapFunction(earliestFinish - latestStart));
                        end
                        interactions(numberOfInteractions + 1) = Interaction(latestStart, earliestFinish, obj.plants(i).index, obj.pollinators(j).index);
                        numberOfInteractions = numberOfInteractions + 1;

                        obj.plants(i).degree = obj.plants(i).degree + 1;
                        obj.pollinators(j).degree = obj.pollinators(j).degree + 1;
                    end
                end
            end
            
            obj.interactions = interactions(1:numberOfInteractions);
        end

        function obj = setExtinctionConditions(obj)
            for i = 1:obj.numberOfPollinators
                obj.pollinators(i).extinctionConditions = getExtinctionConditions(obj, ...
                    obj.pollinators(i), ...
                    obj.interactions([obj.interactions.pollinatorIndex] == obj.pollinators(i).index), obj.minimumReproductionTime);

                obj.pollinators(i).minimalExtinctionConditions = getMinimalExtinctionConditions(obj, obj.pollinators(i).extinctionConditions);
            end
        end

        function extinctionConditions = getExtinctionConditions(~, pollinator, pollinatorInteractions, minimumReproductionTime)

            windowStart = pollinator.startTime;
            windowWidth = (pollinator.endTime - windowStart) * minimumReproductionTime;
            windowEnd = windowStart + windowWidth;

            maxWindows = 2 * length(pollinatorInteractions) + 2;
            maxExtinctionConditions = 2 * length(pollinatorInteractions) + 1;

            extinctionConditions = cell(maxWindows, 1);
            extinctionConditionsPerWindow = cell(maxExtinctionConditions, 1);
            numberOfExtinctionConditions = 0;
            numberOfWindows = 0;
            
            startTimes = [pollinatorInteractions.startTime];
            endTimes = [pollinatorInteractions.endTime];
            plantIndices = [pollinatorInteractions.plantIndex];

            if minimumReproductionTime == 0
                extinctionConditions = {{plantIndices}};
                return;
            end

            while true    
                currentStartTime = windowStart;       
                currentEndTime = min([startTimes(startTimes > currentStartTime), endTimes(endTimes > currentStartTime), windowEnd]);
        
                while currentStartTime < windowEnd
        
                    currentInteractions = and(startTimes <= currentStartTime, endTimes >= currentEndTime);
        
                    extinctionConditionsPerWindow{numberOfExtinctionConditions + 1, 1} = plantIndices(currentInteractions);
                    numberOfExtinctionConditions = numberOfExtinctionConditions + 1;
        
                    currentStartTime = currentEndTime;
                    currentEndTime = min([endTimes(endTimes > currentEndTime), startTimes(startTimes > currentEndTime), windowEnd]);
                end
                numberOfWindows = numberOfWindows + 1;
                extinctionConditionsPerWindow = extinctionConditionsPerWindow(1:numberOfExtinctionConditions, 1);
                extinctionConditions{numberOfWindows} = extinctionConditionsPerWindow;
                numberOfExtinctionConditions = 0;

                if windowEnd >= pollinator.endTime
                    break;
                end

                windowEnd = min([endTimes(endTimes > windowEnd), startTimes(startTimes > windowEnd), pollinator.endTime]);
                windowStart = windowEnd - windowWidth;
            end
            extinctionConditions = extinctionConditions(1:numberOfWindows);
        end
        
        function minimalExtinctionConditions = getMinimalExtinctionConditions(~, extinctionConditions)

            minimalExtinctionConditions = cell(length(extinctionConditions), 1);

            for i = 1:length(extinctionConditions)
                windowExtinctionConditions = extinctionConditions{i};
                      
                isConditionMinimal = true(length(windowExtinctionConditions), 1);
                
                for j = 1:length(windowExtinctionConditions)
                    for k = 1:length(windowExtinctionConditions)
                        if j == k
                            continue;
                        end
                        if length(windowExtinctionConditions{j}) > length(windowExtinctionConditions{k}) && all(ismember(windowExtinctionConditions{k}, windowExtinctionConditions{j}))
                            isConditionMinimal(j) = false;
                            break;
                        end
                    end
                    if ~isConditionMinimal(j)
                        continue;
                    end
                end

                minimalExtinctionConditions{i} = windowExtinctionConditions(isConditionMinimal);
            end
            
            conditions = minimalExtinctionConditions{1};

            for i=2:length(minimalExtinctionConditions)
                minimalWindowConditions = minimalExtinctionConditions{i};
                newConditions = cell(length(conditions), length(minimalWindowConditions));
                for j=1:length(minimalWindowConditions)
                    for k=1:length(conditions)
                        newConditions{k, j} = unique([conditions{k}, minimalWindowConditions{j}]);
                    end
                end
                conditions = reshape(newConditions, 1, []);
                [~, sortedIndices] = sort(cellfun(@length, conditions), 'descend');
                conditions = conditions(sortedIndices);
                isMinimal = true(length(conditions), 1);
                for k = 1:length(conditions)
                    for j=k+1:length(conditions)
                        if all(ismember(conditions{j}, conditions{k}))
                            isMinimal(k) = false;
                            break;
                        end
                    end
                end
                conditions = conditions(isMinimal);
            end
                  
            minimalExtinctionConditions = conditions;
        end

        function isExtinct = isPollinatorExtinct(~, plantExtinctions, pollinator)
            arguments
                ~ 
                plantExtinctions
                pollinator Pollinator
            end
            isExtinct = false;

            for i = 1:length(pollinator.minimalExtinctionConditions)
                if all(ismember(pollinator.minimalExtinctionConditions{i}, plantExtinctions))
                    isExtinct = true;
                    return;
                end
            end
        end
    end
end

