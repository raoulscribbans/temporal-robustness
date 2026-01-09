% produces all figures found in the supplementary material. Requires
% external data files for greenland and euppollnet data, and compares those
% temporal networks with simulation. 

%change file paths as appropriate
greenlandDataFilePath = "data/greenlandData.xlsx";
euppollnetDataFilePath = "data/Interaction_data.csv";

euppollnetNetworkIds = ["Isen_Bjerg_2004" "Isen_Bjerg_2005" "Horbylunde_2005" "Skov_Olesen_2004"];
greendlandNetworkIds = ["Greenland_1996" "Greenland_1997"];

networkIds = [greendlandNetworkIds euppollnetNetworkIds];

networks = [getGreenlandNetworks(greenlandDataFilePath) getEuppollnetNetworks(euppollnetDataFilePath, euppollnetNetworkIds)];

%set number of workers to 0 for non-parallel
workers = length(networks);

networkIterations = 100;
pollinatorWidths = [0.001 0.1:0.1:1];

shuffleBreedingFractions = 0:0.1:1;
reducedBreedingFractions = 0:0.5:1;

results = struct("networkId", num2cell(networkIds));

parfor (i = 1:length(results), workers)
    network = networks(i);

    interactionProbability = sum(network.getBinaryAdjacencyMatrix, 'all')/(network.numberOfPlants * network.numberOfPollinators);

    results(i).shuffledRobustness = getShuffledRobustness(network, networkIterations, shuffleBreedingFractions);
    results(i).shuffledRobustnessConstDegree = getShuffledRobustnessConstDegree(network, networkIterations, shuffleBreedingFractions);
    results(i).uniformShuffledRobustness = getUniformShuffledRobustness(network, networkIterations, shuffleBreedingFractions, interactionProbability);

    uniformPlantWidthDistribution = makedist("Normal", mu = mean([network.plants.intervalWidth]), sigma = 0);
    uniformOverlapFunction = @(x)(x>0)*interactionProbability;
    [~, results(i).uniformDegrees, results(i).uniformTimeAveragedDegrees, results(i).unformIndividualRobustnesses] ...
        = getUniformNetworkData(network, networkIterations, reducedBreedingFractions, pollinatorWidths, uniformPlantWidthDistribution, uniformOverlapFunction);

    results(i).inferredPlantWidthDistribution = fitdist(reshape([network.plants.intervalWidth], network.numberOfPlants, 1), 'Beta');

    interactionGradient = calculateInteractionConstant(network);
    results(i).interactionGradient = interactionGradient;
    linearOverlapFunction = @(x)(x * interactionGradient);

    [~, results(i).linearDegrees, results(i).linearTimeAveragedDegrees, results(i).linearIndividualRobustnesses] ...
        = getUniformNetworkData(network, networkIterations, reducedBreedingFractions, pollinatorWidths, results(i).inferredPlantWidthDistribution, linearOverlapFunction);

    results(i).network = network;
end

skovNetwork = networks(networkIds == "Skov_Olesen_2004");

interactionProbabilityPlot(skovNetwork);

widthDistributionPlot(skovNetwork);

for i = 1:length(results)

    robustness = calculateRobustness(results(i).network, shuffleBreedingFractions);

    formattedTitle = replace(networkIds(i), "_", "\_");

    shuffleFig = shufflePlot( ...
        robustness, ...
        results(i).shuffledRobustness, ...
        results(i).shuffledRobustnessConstDegree, ...
        results(i).uniformShuffledRobustness, ...
        shuffleBreedingFractions);

    title(formattedTitle);

    linearFig = linearInteractionPlot( ...
        results(i).network, ...
        pollinatorWidths, ...
        results(i).unformIndividualRobustnesses, ...
        results(i).uniformDegrees, ...
        results(i).uniformTimeAveragedDegrees, ...
        results(i).linearIndividualRobustnesses, ...
        results(i).linearDegrees, ...
        results(i).linearTimeAveragedDegrees);

    sgtitle(formattedTitle);
end

function networks = getGreenlandNetworks(dataFilePath)
    plantData = readtable(dataFilePath, Sheet = "plants");
    pollinatorDataYear1 = readtable(dataFilePath, Sheet = "pollinators year 1");
    pollinatorDataYear2 = readtable(dataFilePath, Sheet = "pollinators year 2");
    interactionData = readtable(dataFilePath, Sheet = "interactions");
    
    
    %start and end dates are normalised that the season begins on day 1 in data
    %season ends on day ~47 (year 1), and day ~69 (year 2)
    plantNames = string(plantData.Var1);
    pollinatorNamesYear1 = string(pollinatorDataYear1.standardised_name);
    pollinatorNamesYear2 = string(pollinatorDataYear2.standardised_name);
    
    year1Length = max([str2double(plantData.xEnd); str2double(pollinatorDataYear1.xEnd)]);
    year2Length = max([str2double(plantData.xEnd_1); str2double(pollinatorDataYear2.xEnd)]);
    
    buffer = 10^-6;
    
    for i = length(plantNames):-1:1
        plantsYear1(i) = Species();
    
        startTime = (str2double(plantData(i,:).start) - 1) / year1Length;
        endTime = (str2double(plantData(i,:).xEnd)) / year1Length;
    
        plantsYear1(i).time = (endTime + startTime) / 2;
        plantsYear1(i).intervalWidth = (endTime - startTime) + buffer;
    
        plantsYear1(i).name = plantNames(i);
    
        plantsYear2(i) = Species();
    
        startTime = (str2double(plantData(i,:).start_1) - 1) / year2Length;
        endTime = (str2double(plantData(i,:).xEnd_1))/ year2Length;
        plantsYear2(i).time = (endTime + startTime) / 2;
        plantsYear2(i).intervalWidth = (endTime - startTime) + buffer;
    
        plantsYear2(i).name = plantNames(i);
    end
    
    for i = length(pollinatorNamesYear1):-1:1
        pollinatorsYear1(i) = Pollinator();
    
        startTime = (str2double(pollinatorDataYear1(i,:).start) - 1) / year1Length;
        endTime = (str2double(pollinatorDataYear1(i,:).xEnd)) / year1Length;
    
        pollinatorsYear1(i).time = (endTime + startTime) / 2;
        pollinatorsYear1(i).intervalWidth = (endTime - startTime);
    
        pollinatorsYear1(i).name = pollinatorNamesYear1(i);
    end
    
    for i = length(pollinatorNamesYear2):-1:1
        pollinatorsYear2(i) = Pollinator();
    
        startTime = (pollinatorDataYear2(i,:).start - 1) / year2Length;
        endTime = pollinatorDataYear2(i,:).xEnd / year2Length;
    
        pollinatorsYear2(i).time = (endTime + startTime) / 2;
        pollinatorsYear2(i).intervalWidth = (endTime - startTime);
    
        pollinatorsYear2(i).name = pollinatorNamesYear2(i);
    end
    
    adjacencyMatrixYear1 = false(length(plantNames), length(pollinatorNamesYear1));
    adjacencyMatrixYear2 = false(length(plantNames), length(pollinatorNamesYear2));
    
    for i = 1:length(plantNames)
        for j = 1:length(pollinatorNamesYear1)
            interacts = logical(table2array(interactionData(interactionData.Var1 == pollinatorNamesYear1(j), regexprep(plantNames(i),"\s([a-z])", "${upper($1)}"))));
            if any(interacts)
                adjacencyMatrixYear1(i,j) = interacts;
            else
                adjacencyMatrixYear1(i,j) = false;
            end
        end
         for j = 1:length(pollinatorNamesYear2)
            interacts = logical(table2array(interactionData(interactionData.Var1 == pollinatorNamesYear2(j), regexprep(plantNames(i),"\s([a-z])", "${upper($1)}"))));
            if any(interacts)
                adjacencyMatrixYear2(i,j) = interacts;
            else
                adjacencyMatrixYear2(i,j) = false;
            end
        end
    end
    
    networkYear1 = TemporalNetwork(plantsYear1, pollinatorsYear1, 0, 0).setInteractions(adjacencyMatrixYear1);
    networkYear2 = TemporalNetwork(plantsYear2, pollinatorsYear2, 0, 0).setInteractions(adjacencyMatrixYear2);

    networks = [networkYear1, networkYear2];
end

function networks = getEuppollnetNetworks(dataFilePath, networkIds)

    data = readtable(dataFilePath);
    data = data(ismember(data.Network_id, networkIds), :);
    networkSummaryData = getEuppollnetSummaryData(data);

    for i = 1:length(networkIds)
    
        networkId = networkIds(i);
    
        networkData = data(data.Network_id == networkId, :);
        summaryData = networkSummaryData(networkSummaryData.networkId == networkId, :);
    
        seasonLength = days(summaryData.endDate - summaryData.startDate) + 1;
    
        plantNames = string(unique(networkData.Plant_accepted_name));
        pollinatorNames = string(unique(networkData.Pollinator_accepted_name));
    
        plants = repmat(Species(), summaryData.numberOfPlants, 1);
        pollinators = repmat(Pollinator(), summaryData.numberOfPollinators, 1);
    
        adjacencyMatrix = false(summaryData.numberOfPlants, summaryData.numberOfPollinators);
    
        for j = 1:length(plants)
            plantInteractions = networkData(networkData.Plant_accepted_name == plantNames(j), :);
            plantStartDates(j, 1) = min(plantInteractions.Date);
            plantEndDates(j, 1) = max(plantInteractions.Date);
    
            intervalWidth = (days(plantEndDates(j, 1) - plantStartDates(j, 1)) + 1) / seasonLength;
            intervalMidpoint = (days(plantStartDates(j, 1) - summaryData.startDate)) / seasonLength + 0.5 * intervalWidth;
    
            plants(j).intervalWidth = intervalWidth;
            plants(j).time = intervalMidpoint;
    
            for k = 1:length(pollinators)
                if any(plantInteractions.Pollinator_accepted_name == pollinatorNames(k))
                    adjacencyMatrix(j, k) = true;
                end
            end
        end
    
        for j = 1:length(pollinators)
            pollinatorInteractions = networkData(networkData.Pollinator_accepted_name == pollinatorNames(j), :);
            pollinatorStartDates(j, 1) = min(pollinatorInteractions.Date);
            pollinatorEndDates(j, 1) = max(pollinatorInteractions.Date);
    
            intervalWidth = (days(pollinatorEndDates(j, 1) - pollinatorStartDates(j, 1)) + 1) / seasonLength;
            intervalMidpoint = (days(pollinatorStartDates(j, 1) - summaryData.startDate)) / seasonLength + 0.5 * intervalWidth;
    
            pollinators(j).intervalWidth = intervalWidth;
            pollinators(j).time = intervalMidpoint;
    
            pollinators(j).order = string(pollinatorInteractions(1,:).Pollinator_order);
            pollinators(j).genus = string(pollinatorInteractions(1,:).Pollinator_genus);
            pollinators(j).family = string(pollinatorInteractions(1,:).Pollinator_family);
            pollinators(j).name = string(pollinatorInteractions(1,:).Pollinator_accepted_name);
        end

        networks(i) = TemporalNetwork(plants, pollinators, 0, 0).setInteractions(adjacencyMatrix);
    end
end

function networkSummaryData = getEuppollnetSummaryData(data)
    networkUniqueId = unique(string(data.Study_id) + ": " + string(data.Network_id) + " y=" + year(data.Date));
    numberOfNetworks = length(networkUniqueId);
    numberOfDates = zeros(numberOfNetworks, 1);
    numberOfPlants = zeros(numberOfNetworks, 1);
    numberOfPollinators = zeros(numberOfNetworks, 1);
    numberOfInteractions = zeros(numberOfNetworks, 1);
    numberOfUniqueInteractions = zeros(numberOfNetworks, 1);
    
    for i = numberOfNetworks:-1:1
        networkData = data(string(data.Study_id) + ": " + string(data.Network_id) + " y=" + year(data.Date) == networkUniqueId(i), :);
        networkId(i, 1) = networkData(1,:).Network_id;
        studyId(i, 1) = networkData(1,:).Study_id;
        studyYear(i, 1) = year(networkData(1,:).Date);
        numberOfDates(i) = length(unique(networkData.Date));
        startDate(i, 1) = min(networkData.Date);
        endDate(i, 1) = max(networkData.Date);
        numberOfPlants(i) = length(unique(networkData.Plant_accepted_name));
        numberOfPollinators(i) = length(unique(networkData.Pollinator_accepted_name));
        numberOfInteractions(i) = height(networkData);
        numberOfUniqueInteractions(i) = length(unique(string(networkData.Plant_accepted_name) + "+" + string(networkData.Pollinator_accepted_name))); % use '+' as delimitter
    end
    
    networkSummaryData = table(networkUniqueId, networkId, studyId, studyYear, numberOfDates, startDate, endDate, numberOfPlants, numberOfPollinators, numberOfInteractions, numberOfUniqueInteractions);
end

function shuffledRobustness = getShuffledRobustness(network, iterations, breedingFractions)

    shuffledRobustness = zeros(iterations, length(breedingFractions));
    
    overlapIndices = find(network.getOverlapMatrix());
    numberOfInteractions = sum(network.getBinaryAdjacencyMatrix(), 'all');
    
    for j = 1:length(breedingFractions)
        breedingFraction = breedingFractions(j);
        for i = 1:iterations
            interactionIndices = overlapIndices(randperm(length(overlapIndices), numberOfInteractions));
            randomAdjacencyMatrix = false(network.numberOfPlants, network.numberOfPollinators);
            randomAdjacencyMatrix(interactionIndices) = true;
            shuffledRobustness(i, j) = TemporalNetwork([network.plants], [network.pollinators], 0, breedingFraction).setInteractions(randomAdjacencyMatrix).calculateRobustness();
        end
    end
end

function shuffledRobustness = getShuffledRobustnessConstDegree(network, iterations, breedingFractions)

    shuffledRobustness = zeros(iterations, length(breedingFractions));

    adjacencyMatrix = network.getBinaryAdjacencyMatrix();
    overlapMatrix = network.getOverlapMatrix();
    
    for j = 1:length(breedingFractions)
        breedingFraction = breedingFractions(j);
        for i = 1:iterations
            randomAdjacencyMatrix = false(network.numberOfPlants, network.numberOfPollinators);
            for k = 1:network.numberOfPollinators
                pollinatorInteractions = sum(adjacencyMatrix(:, k));
                pollinatorOverlapIndices = find(overlapMatrix(:, k));
    
                randomInteractionIndices = pollinatorOverlapIndices(randperm(length(pollinatorOverlapIndices), pollinatorInteractions));
                randomAdjacencyMatrix(randomInteractionIndices, k) = true;
            end
            shuffledRobustness(i, j) = TemporalNetwork([network.plants], [network.pollinators], 0, breedingFraction).setInteractions(randomAdjacencyMatrix).calculateRobustness();
        end
    end
end
    
function uniformRobustness = getUniformShuffledRobustness(network, iterations, breedingFractions, interactionProbability)
    uniformRobustness = zeros(length(breedingFractions), iterations);
    for j = length(breedingFractions):-1:1
        uniformNetworksShuffle = network.getUniformNullModelNetwork(iterations, mean([network.pollinators.intervalWidth]), breedingFractions(j), true, @(x)((x>0) * interactionProbability));
        for k = 1:iterations
            uniformRobustness(j, k) = uniformNetworksShuffle(k).calculateRobustness();
        end
    end
end

function [uniformNetworks, uniformDegrees, uniformTimeAveragedDegrees, individualRobustnessesUniform] = getUniformNetworkData(network, iterations, breedingFractions, pollinatorWidths, plantWidthDistribution, overlapFunction)
arguments
    network 
    iterations 
    breedingFractions 
    pollinatorWidths 
    plantWidthDistribution = makedist("Normal", mu = mean([network.plants.intervalWidth]), sigma = 0)
    overlapFunction = @(x)(x>0)*(sum(network.getBinaryAdjacencyMatrix, 'all')/(network.numberOfPlants * network.numberOfPollinators))
end
    uniformDegrees = zeros(iterations, length(pollinatorWidths), length(breedingFractions), network.numberOfPollinators);
    uniformTimeAveragedDegrees = zeros(size(uniformDegrees));
    individualRobustnessesUniform = zeros(size(uniformDegrees)); 
    uniformNetworks = repmat(network, length(pollinatorWidths), length(breedingFractions), iterations);

    for j = length(pollinatorWidths):-1:1
        for k = length(breedingFractions):-1:1
            uniformNetworks(j,k,:) = network.getUniformNullModelNetwork( ...
                iterations, ...
                pollinatorWidths(j), ...
                breedingFractions(k), ...
                true, ...
                overlapFunction, ...
                plantWidthDistribution);
            for l = 1:iterations
                uniformDegrees(l,j,k,:) = [uniformNetworks(j,k,l).pollinators.degree];
                if pollinatorWidths(j) < 0.01
                    uniformTimeAveragedDegrees(l,j,k,:) = uniformNetworks(j,k,l).getTimeAveragedDegrees(0.0001);
                else
                    uniformTimeAveragedDegrees(l,j,k,:) = uniformNetworks(j,k,l).getTimeAveragedDegrees(0.001);
                end
                individualRobustnessesUniform(l,j,k,:) = uniformNetworks(j,k,l).calculateIndividualRobustness();
            end
        end
    end
end

function interactionGradient = calculateInteractionConstant(network)
    doesInteract = double(reshape(network.getBinaryAdjacencyMatrix(), network.numberOfPlants * network.numberOfPollinators, 1));
    overlapAmount = reshape(network.getWeightedAdjacencyMatrix(), network.numberOfPlants * network.numberOfPollinators, 1);

    %use linear fit to obtain constant of proportionality
    interactionGradient = plotHelper.linearFit(overlapAmount, doesInteract, false).gradient;
end

function robustness = calculateRobustness(network, breedingFractions)
    robustness = zeros(length(breedingFractions), 1);
    for i = 1:length(breedingFractions)
        robustness(i) = ...
            TemporalNetwork([network.plants], [network.pollinators], 0, breedingFractions(i)) ...
            .setInteractions(network.getBinaryAdjacencyMatrix) ...
            .calculateRobustness();
    end
end

function f = shufflePlot(robustness, shuffledRobustness, shuffledRobustnessConstDegree, uniformShuffledRobustness, breedingFractions)
    f = figure();

    co = colororder;

    p(1) = plot(breedingFractions, robustness, LineWidth=2);
    
    hold on;
    p(2) = plot(breedingFractions, mean(shuffledRobustness, 1), LineWidth=2, LineStyle="--");
    hold on;
    p(3) = plot(breedingFractions, mean(shuffledRobustnessConstDegree, 1), LineWidth=2, LineStyle="-.");
    hold on;
    p(4) = plot(breedingFractions, mean(uniformShuffledRobustness, 2), LineWidth=2, LineStyle=":"); %mean over dim 2 bc i did it wrong before
    
    legend(p, ["data" "shuffled (const. connectance)" "shuffled (const. degrees)" "uniform model"]);
    
    xlabel("Minimum Reproduction Time a");
    ylabel("Robustness R");
        
    xlim([0 1]);
    ylim([0 1]);
end

function f = linearInteractionPlot(network, modelPollinatorWidths, unformIndividualRobustnesses, uniformDegrees, uniformTimeAveragedDegrees, linearIndividualRobustnesses, linearDegrees, linearTimeAveragedDegrees)
    f = figure();
    title("tets")
    co = colororder;

    tl = tiledlayout(3, 2);

    nexttile(1, [2 1])
    
    individualAggregatedRobustnesses = TemporalNetwork([network.plants], [network.pollinators], 0, 0) ...
        .setInteractions(network.getBinaryAdjacencyMatrix()) ...
        .calculateIndividualRobustness();
   
    uniformVsLinearPlot( ...
        [network.pollinators.intervalWidth], ...
        individualAggregatedRobustnesses, ...
        modelPollinatorWidths, ...
        mean(unformIndividualRobustnesses(:,:,1,:), [1, 4]), ...
        sqrt(var(unformIndividualRobustnesses(:,:,1,:), 0, [1, 4])), ...
        mean(linearIndividualRobustnesses(:,:,1,:), [1, 4]), ...
        sqrt(var(linearIndividualRobustnesses(:,:,1,:), 0, [1, 4])), ...
        true);
    
    title("Time Aggregated (a=0)");

    nexttile([2 1])
    
    
   individualTemporalRobustnesses = TemporalNetwork([network.plants], [network.pollinators], 0, 1) ...
        .setInteractions(network.getBinaryAdjacencyMatrix()) ...
        .calculateIndividualRobustness();

    uniformVsLinearPlot( ...
        [network.pollinators.intervalWidth], ...
        individualTemporalRobustnesses, ...
        modelPollinatorWidths, ...
        mean(unformIndividualRobustnesses(:,:,end,:), [1, 4]), ...
        sqrt(var(unformIndividualRobustnesses(:,:,end,:), 0, [1, 4])), ...
        mean(linearIndividualRobustnesses(:,:,end,:), [1, 4]), ...
        sqrt(var(linearIndividualRobustnesses(:,:,end,:), 0, [1, 4])), ...
        false);

    title("Temporal (a=1)");

    nexttile([1 1]);
    
    degreePlot( ...
        [network.pollinators.intervalWidth], ...
        sum(network.getBinaryAdjacencyMatrix, 1), ...
        modelPollinatorWidths, ...
        mean(uniformDegrees, [1, 3, 4]), ...
        mean(linearDegrees, [1, 3, 4]));

    ylabel("Degree k")
    xlabel("Active Period Duration w_A");

    nexttile([1 1]);
    
    degreePlot( ...
        [network.pollinators.intervalWidth], ...
        network.getTimeAveragedDegrees(0.001), ...
        modelPollinatorWidths, ...
        mean(uniformTimeAveragedDegrees, [1, 3, 4]), ...
        mean(linearTimeAveragedDegrees, [1, 3, 4]));
        
    ylabel("Time-Averaged Degree")
    xlabel("Active Period Duration w_A");
    
    tl.Padding = "compact";
    %tl.TileSpacing = 'compact';
    
    f.Position(3:4)=[3, 1.1] * 500;
    
    fontsize(14, 'points')
end

function p = uniformVsLinearPlot(pollinatorWidths, individualRobustnesses, modelPollinatorWidths, uniformIndividualRobustnessesMean, uniformIndividualRobustnessesSd, linearIndividualRobustnessesMean, linearIndividualRobustnessesSd, showLegend)
        
        co = colororder;
        
        s = scatter(pollinatorWidths, individualRobustnesses, 48, co(1, :), Marker = "square", LineWidth=2);
        hold on;
        
        p(1) = plot(modelPollinatorWidths, uniformIndividualRobustnessesMean, LineWidth=2, Color=co(3,:));
        hold on;
        plotHelper.fillError(modelPollinatorWidths, ...
            uniformIndividualRobustnessesMean - uniformIndividualRobustnessesSd, ...
            uniformIndividualRobustnessesMean + uniformIndividualRobustnessesSd, ...
            co(3, :));

        hold on;
        
        p(2)= plot(modelPollinatorWidths, linearIndividualRobustnessesMean, LineWidth=2, Color=co(4,:));
        hold on;
        plotHelper.fillError(modelPollinatorWidths, ...
            linearIndividualRobustnessesMean - linearIndividualRobustnessesSd, ...
            linearIndividualRobustnessesMean + linearIndividualRobustnessesSd, ...
            co(4, :));
        hold on;
        
        xlim([0 1]);
        ylim([0 1]);
        
        ylabel("Individual Robustness r")
        if showLegend
            legend([s, p(1), p(2)], ["Data", "Uniform Model", "Linear Interaction Model"], location = "southeast");
        end
end

function degreePlot(pollinatorWidths, degrees, modelPollinatorWidths, uniformDegrees, linearDegrees)

        co = colororder;

        scatter(pollinatorWidths, degrees, LineWidth=2);
        hold on;
        
        plot(modelPollinatorWidths, uniformDegrees, Color=co(3, :), LineWidth=2);
        
        hold on;
        plot(modelPollinatorWidths, linearDegrees, Color=co(4, :), LineWidth = 2);
       
end

function interactionProbabilityPlot(network)
    overlapAmouts = reshape(network.getOverlapMatrix(), [], 1);
    interacts = double(reshape(network.getBinaryAdjacencyMatrix(), [], 1));
    
    binStep = 0.05;
    bins = 0:binStep:1;
    
    figure();
    
    scatter(overlapAmouts, interacts, LineWidth=1, Marker="x");
    hold on;

    y = zeros(length(bins) - 1, 1);
    
    for i = 1:length(y)
        overlapsInbin = overlapAmouts >= bins(i) & overlapAmouts < bins(i + 1);
        y(i) = sum(interacts .* overlapsInbin, 'all') / sum(overlapsInbin, 'all');
    end
    
    x = (bins(2:end) + bins(1:end-1)) / 2;
    
    scatter(x, y, LineWidth=2);
    
    hold on;
    
    interactionConstant = plotHelper.linearFit(overlapAmouts, interacts, false).gradient;
    
    plot(0:1, interactionConstant * (0:1), LineWidth=2, LineStyle="--");
    
    ylim([0 1]);
    xlim([0 1]);
    
    legend("Data", "Data (Binned)", "Linear Fit", Location="northwest")
    
    xlabel("Overlap Duration");
    ylabel("Interaction Probability");
    
    fontsize(14, 'points')
end

function widthDistributionPlot(network)
    plantWidths = [network.plants.intervalWidth];
    
    figure();
    
    histogram(plantWidths, 0:0.1:1, Normalization="pdf");
    hold on;
    
    fitDistribution = fitdist(plantWidths', 'beta');
    xValues = linspace(0, 1, 100);
    plot(xValues, fitDistribution.pdf(xValues), LineWidth=2, LineStyle="--");
    
    xlim([0 1]);
    xlabel("Plant Active Period Duration W_P")
    ylabel("Probability Density");
    
    legend("Data", "Beta Dist.", Location="northeast")
    
    fontsize(14, 'points')
end