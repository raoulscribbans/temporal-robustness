% reads data located at https://zenodo.org/records/17123949 into a temporal
% network. Evalutes robustness and compares it to various null models.
% Prodeces figures 2 and 4 of the main text

zenodoDoi = "17123949";
fileName = "2_Overlap_time.csv";

data = getZenodoData(zenodoDoi, fileName);

years = unique(data.Year);

pollinatorNames = unique(string(data.Pollinator_gen_sp));
plantNames = unique(string(data.Plant_gen_sp));

numberOfPlants = length(plantNames);
numberOfPollinators = length(pollinatorNames);

networks = repmat(TemporalNetwork(Species(), Pollinator(), 0, 1), length(years), 1);

breedingFractions = 0:0.1:1;

robustnesses = zeros(length(years), length(breedingFractions));
individualRobustnesses = zeros(length(years), length(breedingFractions), numberOfPollinators);

pollinatorWidths = zeros(length(years), numberOfPollinators);
plantWidths = zeros(length(years), numberOfPlants);

overlaps = false(length(years), numberOfPlants, numberOfPollinators);
overlapAmount = zeros(length(years), numberOfPlants, numberOfPollinators);

interactionProbability = zeros(length(years), 1);

shuffleIterations = 100;
uniformIterations = 100;

shuffledRobustnesses = zeros(length(years), shuffleIterations, length(breedingFractions));
shuffledRobustnessesPol = zeros(length(years), shuffleIterations, length(breedingFractions));
uniformModelRobustnesses = zeros(length(years), shuffleIterations, length(breedingFractions));

%Time-averaged degrees
timeStep = 0.0001;
times = 0:timeStep:1;

timeAveragedDegrees = zeros(length(years), numberOfPollinators);
timeAveragedDegreesSquared = zeros(length(years), numberOfPollinators);

%small buffer to account for rounding errors
%added to each plant interval
buffer = 10^-6;

for i = 1:1:length(years)
    yearlyData = data(data.Year == years(i),:);

    adjacencyMatrix = false(numberOfPlants, numberOfPollinators);
    
    for ii = 1:numberOfPlants
        for j = 1:numberOfPollinators
            adjacencyMatrix(ii, j) = any(yearlyData.Plant_gen_sp == plantNames(ii) & yearlyData.Pollinator_gen_sp == pollinatorNames(j));
        end
    end

    seasonStart = min([yearlyData.First_seenPoll; yearlyData.First_seenPlan]);
    seasonEnd = max([yearlyData.Last_seenPoll; yearlyData.Last_seenPlan]);

    seasonLength = seasonEnd - seasonStart + 1;

    for j = numberOfPollinators:-1:1
        pollinators(j) = Pollinator();

        %take min to ensure single value
        startDay = min(table2array(yearlyData(yearlyData.Pollinator_gen_sp == pollinatorNames(j), "First_seenPoll"))); 
        endDay = min(table2array(yearlyData(yearlyData.Pollinator_gen_sp == pollinatorNames(j), "Last_seenPoll")));
        
        pollinators(j).intervalWidth = (endDay + 1 - startDay) / seasonLength;
        pollinators(j).time = (startDay - seasonStart) / seasonLength + pollinators(j).intervalWidth / 2;
        pollinators(j).name = pollinatorNames(j);

        pollinatorWidths(i, j) = pollinators(j).intervalWidth;
    end

    for j = numberOfPlants:-1:1
        plants(j) = Species();

        %take min to ensure single value
        startDay = min(table2array(yearlyData(yearlyData.Plant_gen_sp == plantNames(j), "First_seenPlan")));
        endDay = min(table2array(yearlyData(yearlyData.Plant_gen_sp == plantNames(j), "Last_seenPlan")));

        plants(j).intervalWidth = (endDay + 1 - startDay) / seasonLength;
        plants(j).time = (startDay - seasonStart) / seasonLength + plants(j).intervalWidth / 2;
        plants(j).intervalWidth = plants(j).intervalWidth + buffer;
        plants(j).name = plantNames(j);

        plantWidths(i, j) = plants(j).intervalWidth;
    end
    
    for j = 1:numberOfPlants
        for k = 1:numberOfPollinators
            overlaps(i, j, k) = plants(j).startTime < pollinators(k).endTime && plants(j).endTime > pollinators(k).startTime;
            if overlaps(i, j, k)
                sortedTimes = sort([plants(j).startTime, plants(j).endTime, pollinators(k).startTime, pollinators(k).endTime]);
                overlapAmounts(i, j, k) = sortedTimes(3) - sortedTimes(2); %overlap is given by difference of two centremost times
            end
        end
    end

    overlapMatrix = reshape(overlaps(i,:,:), size(overlaps, 2), size(overlaps, 3));

    interactionProbability(i) = sum(overlapMatrix & adjacencyMatrix, 'all') / sum(overlapMatrix, 'all');
    
    for j = 1:length(breedingFractions)
        robustnesses(i,j) = TemporalNetwork(plants, pollinators, 0, breedingFractions(j)).setInteractions(adjacencyMatrix).calculateRobustness;
        individualRobustnesses(i,j,:) = TemporalNetwork(plants, pollinators, 0, breedingFractions(j)).setInteractions(adjacencyMatrix).calculateIndividualRobustness;
    end

    overlapIndices = find(overlaps(i, :, :));
    interactionMatrix = overlapMatrix & adjacencyMatrix;
    numberOfinteractions = sum(interactionMatrix, 'all');

    %randomly shuffle interaction, conserve total number of interations
    for j = 1:length(breedingFractions)
        breedingFraction = breedingFractions(j);
        for l = 1:shuffleIterations
            interactionIndices = overlapIndices(randperm(length(overlapIndices), numberOfinteractions));
            randomAdjacencyMatrix = false(numberOfPlants, numberOfPollinators);
            randomAdjacencyMatrix(interactionIndices) = true;
            shuffledRobustnesses(i, l, j) = TemporalNetwork(plants, pollinators, 0, breedingFraction).setInteractions(randomAdjacencyMatrix).calculateRobustness();
        end
    end

    %randomly shuffle interaction, conserve degree of each pollinator (no
    %effect for time aggregated)
    for j = 1:length(breedingFractions)
        breedingFraction = breedingFractions(j);
        for l = 1:shuffleIterations
            randomAdjacencyMatrix = false(numberOfPlants, numberOfPollinators);
            for k = 1:length(pollinators)
                pollinatorInteractions = sum(interactionMatrix(:, k));
                pollinatorOverlapIndices = find(overlapMatrix(:, k));

                randomInteractionIndices = pollinatorOverlapIndices(randperm(length(pollinatorOverlapIndices), pollinatorInteractions));
                randomAdjacencyMatrix(randomInteractionIndices, k) = true;
            end
            shuffledRobustnessesPol(i, l, j) = TemporalNetwork(plants, pollinators, 0, breedingFraction).setInteractions(randomAdjacencyMatrix).calculateRobustness();
        end
    end

    % %uniform network model using mean interval width
    pollinatorWidthDistribution = makedist("Normal", mu = mean([pollinators.intervalWidth]), sigma = 0) ;
    plantWidthDistribution = makedist("Normal", mu = mean([plants.intervalWidth]), sigma = 0);
    % 
    for j = 1:length(breedingFractions)
        for l = 1:shuffleIterations
            breedingFraction = breedingFractions(j);

            for k = length(pollinators):-1:1
                modelPollinators(k) = Pollinator(pollinatorWidthDistribution, true);
            end
            for k = length(plants):-1:1
                modelPlants(k) = Species(plantWidthDistribution, true);
            end

            uniformModelRobustnesses(i, l, j) = TemporalNetwork(modelPlants, modelPollinators, interactionProbability(i), breedingFraction) ...
                .calculateRobustness();


        end
    end

    timedDegrees = zeros(length(times), numberOfPollinators);
    
    for j = 1:length(times)       
        plantsActive = [plants.startTime] < times(j) & [plants.endTime] > times(j);
        pollinatorsActive = [pollinators.startTime] < times(j) & [pollinators.endTime] > times(j);

        timedDegrees(j, :) = (plantsActive * adjacencyMatrix) .* pollinatorsActive;
        timeAveragedDegrees(i, :) = timeAveragedDegrees(i, :) + timedDegrees(j, :);
        timeAveragedDegreesSquared(i, :) = timeAveragedDegreesSquared(i, :) + timedDegrees(j, :).^2;
    end

    timeAveragedDegrees(i, :) = timeAveragedDegrees(i, :) ./ ([pollinators.intervalWidth]) / length(times);
    timeAveragedDegreesSquared(i, :) = timeAveragedDegreesSquared(i, :) ./ ([pollinators.intervalWidth]) / length(times);
end


meanPlantWidth = mean(mean(plantWidths));
uniformPollinatorWidths = [0.000001 0.05:0.05:1]; %small non-zero value to approximate 0 (having w=0 leads to 0 interactions)

uniformTimeAveragedDegrees = zeros(uniformIterations, length(uniformPollinatorWidths));
unformDegrees = zeros(uniformIterations, length(uniformPollinatorWidths));


plantWidthDistribution = makedist("Normal", mu = meanPlantWidth, sigma = 0);
inferredPlantWidthDistribution = fitdist(mean(plantWidths, 1)', "Beta");

breedingFractionsAggOrTemp = [0 1];

meanOverlapAmount = reshape(mean(overlapAmounts, 1), numberOfPlants, numberOfPollinators);
meanInteractionProbability = reshape(mean(overlaps, 1), numberOfPlants, numberOfPollinators) .* adjacencyMatrix;

linearFitCoefficients = plotHelper.linearFit(reshape(meanOverlapAmount, [], 1), reshape(meanInteractionProbability, [], 1), false);

linearOverlapFunction = @(x) x * linearFitCoefficients.gradient;

uniformNetworks = repmat(TemporalNetwork(Species(),Pollinator()), uniformIterations, length(uniformPollinatorWidths), length(breedingFractionsAggOrTemp));
linearInteractionNetworks = repmat(TemporalNetwork(Species(),Pollinator()), uniformIterations, length(uniformPollinatorWidths), length(breedingFractionsAggOrTemp));

for i = 1:length(breedingFractionsAggOrTemp)
    for j = 1:length(uniformPollinatorWidths)
        pollinatorWidthDistribution = makedist("Normal", mu = uniformPollinatorWidths(j), sigma = 0);
        relativeOverlapFunction = @(x) x * linearFitCoefficientsRelativeOverlap.gradient / uniformPollinatorWidths(k);

        uniformNetworks(:, j, i) = TemporalNetwork.createFromDistributions( ...
            plantWidthDistribution, ...
            pollinatorWidthDistribution, ...
            numberOfPlants, ...
            numberOfPollinators, ...
            true, ...
            [uniformIterations, 1], ...
            mean(interactionProbability), ...
            breedingFractionsAggOrTemp(i));

        linearInteractionNetworks(:, j, i) = TemporalNetwork.createFromDistributions( ...
            inferredPlantWidthDistribution, ...
            pollinatorWidthDistribution, ...
            numberOfPlants, ...
            numberOfPollinators, ...
            true, ...
            [uniformIterations, 1], ...
            0, ...
            breedingFractionsAggOrTemp(i), ...
            linearOverlapFunction);
    end
end

uniformIndividualRobustness = zeros([size(uniformNetworks), length(pollinators)]);
uniformDegrees = uniformIndividualRobustness;
uniformTimeAvDegrees = uniformIndividualRobustness;

for i = 1:size(uniformIndividualRobustness, 1)
    for j = 1:size(uniformIndividualRobustness, 2)
        for k = 1:size(uniformIndividualRobustness, 3)
            uniformIndividualRobustness(i,j,k,:) = uniformNetworks(i,j,k).calculateIndividualRobustness();
            uniformDegrees(i,j,k,:) = [uniformNetworks(i,j,k).pollinators.degree];
            if j == 1
                timeStep = 0.00001;
            else 
                timeStep = 0.001;
            end
            uniformTimeAvDegrees(i,j,k,:) = uniformNetworks(i,j,k).getTimeAveragedDegrees(timeStep);
        end
    end
end

uniformRobustness = mean(uniformIndividualRobustness, 4);

linearInteractionIndividualRobustness = zeros([size(uniformNetworks), length(pollinators)]);
linearInteractionDegrees = uniformIndividualRobustness;
linearInteractionTimeAvDegrees = uniformIndividualRobustness;

for i = 1:size(linearInteractionIndividualRobustness, 1)
    for j = 1:size(linearInteractionIndividualRobustness, 2)
        for k = 1:size(linearInteractionIndividualRobustness, 3)
            linearInteractionIndividualRobustness(i,j,k,:) = linearInteractionNetworks(i,j,k).calculateIndividualRobustness();
            linearInteractionDegrees(i,j,k,:) = [linearInteractionNetworks(i,j,k).pollinators.degree];
            if j == 1
                timeStep = 0.00001;
            else 
                timeStep = 0.001;
            end
            linearInteractionTimeAvDegrees(i,j,k,:) = linearInteractionNetworks(i,j,k).getTimeAveragedDegrees(timeStep);
        end
    end
end

linearInteractionRobustness = mean(linearInteractionIndividualRobustness, 4);


aggregatedRobustness = squeeze(mean(individualRobustnesses(:, 1, :), 1));
temporalRobustness = squeeze(mean(individualRobustnesses(:, end, :), 1));
temporalRobustnessAHalf = squeeze(mean(individualRobustnesses(:, breedingFractions == 0.5, :), 1));
meanWidths = mean(pollinatorWidths, 1);

f = figure();

tl = tiledlayout(6, 3);

%time plot
nexttile([6 1]);
    
    co = colororder;
    
    [~, orderPlantIndices] = sort([plants.time], 'ascend');
    [~, orderPollinatorsIndices] = sort([pollinators.time], 'ascend');
    
    timePlot = plotHelper.timePlot(plants(orderPlantIndices), pollinators(orderPollinatorsIndices), 0.75, co(3,:), co(2,:));
    
    xlim([0 1]);
    xlabel("Active Period");

    t(1) = title('A');

%network
nexttile([3 1]);

    co = colororder;
    
    [orderedPlantTimes, orderPlantIndices] = sort([plants.time], 'descend');
    [orderedPollinatorTimes, orderPollinatorsIndices] = sort([pollinators.time], 'ascend');
    
    weightedAdjacency = reshape(overlapAmounts(end,:,:), numberOfPlants, numberOfPollinators) .* adjacencyMatrix;
    
    bigAdjacencyMatrix = [zeros(numberOfPlants), weightedAdjacency(orderPlantIndices, orderPollinatorsIndices); weightedAdjacency(orderPlantIndices, orderPollinatorsIndices)', zeros(numberOfPollinators)];
    
    g = graph(bigAdjacencyMatrix, [plantNames(orderPlantIndices); pollinatorNames(orderPollinatorsIndices)]);
    
    yData = [ones(numberOfPlants, 1); ones(numberOfPollinators, 1) * 2];
    xData = [1:numberOfPlants, (1:numberOfPollinators) * numberOfPlants / numberOfPollinators];
    
    markerScaleFactor = 10;
    edgeScaleFactor = 5;
    
    h = plot(g, ...,
        "Layout","circle",...
        "NodeColor", [repmat(co(3,:), numberOfPlants, 1); repmat(co(2,:), numberOfPollinators, 1)], ...
        'LineWidth', g.Edges.Weight * edgeScaleFactor, ... 
        "MarkerSize", [plants(orderPlantIndices).intervalWidth, pollinators(orderPollinatorsIndices).intervalWidth] * markerScaleFactor, ...
        Marker=[repmat("o", 1, numberOfPlants), repmat("square", 1, numberOfPollinators)], ...
        NodeLabel = {}, ...
        NodeFontSize=16);

    t(2) = title('B');

%Robustness vs pollinator width
nexttile([3 1])

    for i = 1:length(aggregatedRobustness)
    
        for j = 1:3
            switch j
                case 1
                    s(i, j) = scatter([meanWidths(i)], aggregatedRobustness(i), 48, co(1, :), Marker = "square", LineWidth=2);
                case 2
                    s(i, j) = scatter([meanWidths(i)], temporalRobustnessAHalf(i), 48, co(2, :), Marker = "x", LineWidth=2);
                otherwise
                    s(i, j) = scatter([meanWidths(i)], temporalRobustness(i), 48, co(3, :), Marker = "diamond", LineWidth=2);
            end
            hold on;
        end
        plot([meanWidths(i) meanWidths(i) meanWidths(i)], [aggregatedRobustness(i) temporalRobustnessAHalf(i) temporalRobustness(i)], Color=[0 0 0], LineStyle="--");
        hold on;
    end
    
    legend([s(1, 1) s(1, 2) s(1, 3)], "Aggregated Robustness (a=0)", "Temporal Robustness (a=1/2)", "Temporal Robustness (a=1)", Location = "southwest");
        
    
    xlim([0 1]);
    ylim([0 1]);
    
    ylabel("Individual Robustness r")
    t(3) = title('C');

%shuffle interactions
nexttile([3 1]);
    
    clear p;
    co = colororder;
    p(1) = plot(breedingFractions, mean(robustnesses, 1), LineWidth=2);
    hold on;
    plotHelper.fillError(breedingFractions, min(robustnesses,[],1), max(robustnesses,[],1), co(1,:));
    
    hold on;
    p(2) = plot(breedingFractions, squeeze(mean(mean(shuffledRobustnesses, 1), 2)), LineWidth=2, LineStyle="--");
    hold on;
    p(3) = plot(breedingFractions, squeeze(mean(mean(shuffledRobustnessesPol, 1), 2)), LineWidth=2, LineStyle="-.");
    hold on;
    p(4) = plot(breedingFractions, squeeze(mean(mean(uniformModelRobustnesses, 1), 2)), LineWidth=2, LineStyle=":");
    
    legend(p, ["data" "shuffled (const. connectance)" "shuffled (const. degrees)" "uniform model"], Location = "southwest");
    
    xlabel("Minimum Reproduction Time a");
    ylabel("Robustness R");
    
    xlim([0 1]);
    ylim([0.4 1]);

    t(4) = title('D');

%change in r
nexttile([3 1])

    Y = aggregatedRobustness - temporalRobustness;
    scatter(meanWidths, Y, LineWidth=2);
    hold on;
    X = [ones(length(meanWidths),1) meanWidths'];
    b = X\Y;
    plot([0 1],  [1 0; 1 1]*b, LineStyle="--", LineWidth=2)
    ylim([0 0.35]);
    xlim([0 1]);
    
    xlabel("Active Period Duration w_A");
    ylabel("Change in Individual Robustness \Delta{r}")

    t(5) = title('E');

tl.Padding = "compact";
tl.TileSpacing = 'compact';

f.Position(1:2) = [0, 0];
f.Position(3:4)=[1000, 500] * 1.7;

fontsize(14, 'points');

for i = 1:length(t)
    t(i).FontSize = 30;
end

degrees = zeros(1, numberOfPollinators);

for i = 1:length(years)
    interactions = adjacencyMatrix & reshape(overlaps(i, :, :), numberOfPlants, numberOfPollinators);
    degrees = degrees + sum(interactions) / length(years);
end


f = figure();
tl = tiledlayout(3, 2);
nexttile(1, [2 1])

    aggs = squeeze(mean(individualRobustnesses(:, 1, :), 1));
    temps = squeeze(mean(individualRobustnesses(:, end, :), 1));
    halfs = squeeze(mean(individualRobustnesses(:, breedingFractions == 0.5, :), 1));
    
    co = colororder;
    
    s = scatter(mean(pollinatorWidths, 1), aggs, 48, co(1, :), Marker = "square", LineWidth=2);
    hold on;
    
    p(1) = plot(uniformPollinatorWidths, mean(uniformRobustness(:,:,1), 1), LineWidth=2, Color=co(3,:));
    hold on;
    plotHelper.fillError(uniformPollinatorWidths, ...
        mean(uniformRobustness(:,:,1), 1) - sqrt(var(uniformIndividualRobustness(:,:,1, :), 0,[1, 3, 4])), ...
        mean(uniformRobustness(:,:,1), 1) + sqrt(var(uniformIndividualRobustness(:,:,1, :), 0,[1, 3, 4])), ...
        co(3, :));
    hold on;
    
    p(2) = plot(uniformPollinatorWidths, mean(linearInteractionIndividualRobustness(:,:,1, :), [1, 4]), LineWidth=2, Color=co(4,:));
    hold on;
    plotHelper.fillError(uniformPollinatorWidths, ...
        mean(linearInteractionIndividualRobustness(:,:,1, :), [1, 3, 4]) - sqrt(var(linearInteractionIndividualRobustness(:,:,1, :), 0,[1, 3, 4])), ...
        mean(linearInteractionIndividualRobustness(:,:,1, :), [1, 3, 4]) + sqrt(var(linearInteractionIndividualRobustness(:,:,1, :), 0,[1, 3, 4])), ...
        co(4, :));
    hold on;
    
    xlim([0 1]);
    ylim([0 1]);
    
    ylabel("Individual Robustness r")
    
    legend([s, p(1), p(2)], ["Data", "Uniform Model", "Linear Interaction Model"], location = "southeast");

    title("Time-Aggregated Robustness (a=0)");

nexttile([2 1])

    temps = squeeze(mean(individualRobustnesses(:, end, :), 1));
    s = scatter(mean(pollinatorWidths, 1), temps, 48, co(1, :), Marker = "square", LineWidth=2);
    hold on;
    
    p(1) = plot(uniformPollinatorWidths, mean(uniformRobustness(:,:,2), 1), LineWidth=2, Color=co(3,:));
    hold on;
    plotHelper.fillError(uniformPollinatorWidths, ...
        mean(uniformRobustness(:,:,2), 1) - sqrt(var(uniformIndividualRobustness(:,:,2, :), 0,[1, 4])), ...
        mean(uniformRobustness(:,:,2), 1) + sqrt(var(uniformIndividualRobustness(:,:,2, :), 0,[1, 4])), ...
        co(3, :));
    hold on;
    
    p(2) = plot(uniformPollinatorWidths, mean(linearInteractionIndividualRobustness(:,:,2, :), [1, 4]), LineWidth=2, Color=co(4,:));
    hold on;
    plotHelper.fillError(uniformPollinatorWidths, ...
        mean(linearInteractionIndividualRobustness(:,:,2,:), [1, 4]) - sqrt(var(linearInteractionIndividualRobustness(:,:,2, :), 0,[1, 4])), ...
        mean(linearInteractionIndividualRobustness(:,:,2,:), [1, 4]) + sqrt(var(linearInteractionIndividualRobustness(:,:,2, :), 0,[1, 4])), ...
        co(4, :));
    hold on;
    
    xlim([0 1]);
    ylim([0 1]);
    
    ylabel("Individual Robustness r")
        
    title("Temporal Robustness (a=1)");

nexttile([1 1]);

    s(1) = scatter(mean(pollinatorWidths), degrees, LineWidth=2);
    hold on;
    
    hold on;
    plot(uniformPollinatorWidths, mean(uniformDegrees, [1, 3, 4]), Color=co(3, :), LineWidth=2);
    
    hold on;
    plot(uniformPollinatorWidths, mean(linearInteractionDegrees, [1,3,4]), Color=co(4, :), LineWidth = 2);
    
    ylabel("Degree k")
    xlabel("Active Period Duration w_A");
    
nexttile();

    scatter(mean(pollinatorWidths), mean(timeAveragedDegrees), LineWidth=2);
    
    hold on;
    plot(uniformPollinatorWidths, mean(uniformTimeAvDegrees, [1, 3, 4]), Color=co(3, :), LineWidth = 2);
    
    hold on;
    plot(uniformPollinatorWidths, mean(linearInteractionTimeAvDegrees, [1,3,4]), Color=co(4, :), LineWidth = 2);
    
    ylabel("Time-Averaged Degree")
    xlabel("Active Period Duration w_A");

tl.Padding = "compact";

f.Position(3:4)=[3, 1.1] * 500;

fontsize(14, 'points')

function data = getZenodoData(doi, targetFilename)
% DOWNLOAD_ZENODO_FILE Download a file from Zenodo using its DOI
    fileUrl = "https://zenodo.org/api/records/" + doi + "/files/" + targetFilename + "/content";
    data = readtable(fileUrl);
end