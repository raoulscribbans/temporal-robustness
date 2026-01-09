% Simulates uniform model temporal networks and evaluates their robusntess. 
% Results are plotted to produce figure 3 of the main text. 

numberOfPlants = 200;
numberOfPollinators  = 200;
numberOfProcesses = 40;
plantWidth = 0.05;
interactionProbability = 0.5;
iterations = 100;
extinctionSequences = 100;

%set number of parallel processes - set to 0 for non-parallel
numberOfworkers = 12;

%set variables for r vs theta plot
pollinatorWidths = plantWidth * (0.1:0.1:4);
minimumReproductionTimes = 0:0.2:1; 

dataRobustnessVsTheta = getRobustnessData( ...
    numberOfPlants, ...
    numberOfPollinators, ...
    plantWidth, ...
    interactionProbability, ...
    pollinatorWidths, ...
    minimumReproductionTimes, ...
    iterations, ...
    extinctionSequences, ...
    numberOfworkers);

%set variables for R vs a plot
pollinatorWidths = plantWidth * ([0.5 1 2 3]);
minimumReproductionTimes = 0:0.05:1;

dataRobustnessVsMrt = getRobustnessData( ...
    numberOfPlants, ...
    numberOfPollinators, ...
    plantWidth, ...
    interactionProbability, ...
    pollinatorWidths, ...
    minimumReproductionTimes, ...
    iterations, ...
    extinctionSequences, ...
    numberOfworkers);

f = figure();
tl = tiledlayout(2, 4);

co = colororder;
markers = ["o", "x", "+", "square", "diamond", "^"];

exampleNetwork = TemporalNetwork.createFromDistributions( ...
    makedist('Normal',mu=0.15, sigma=0), ...
    makedist('Normal',mu=0.3, sigma=0), ...
    15, ...
    15, ...
    true, ... 
    1, ...
    0.5);

nexttile([2 1])
plotExampleTime(exampleNetwork);
t(1) = title("A");

nexttile([1 1])
plotExampleNetwork(exampleNetwork)
t(2) = title("B");

nexttile([2 2])   
plotRobustnessVsTheta(dataRobustnessVsTheta, markers);
t(3) = title("C");

nexttile([1 1])
plotRobustnessVsMrt(dataRobustnessVsMrt, markers);
t(4) = title("D");

tl.Padding = "compact";
tl.TileSpacing = 'compact';

f.Position(3:4)=[4, 2] * 400;

fontsize(14, 'points');

for i = 1:length(t)
    t(i).FontSize = 30;
end

function plotExampleTime(network)
    co = colororder();

    plants = [network.plants];
    pollinators = [network.pollinators];
       
    [~, orderPlantIndices] = sort([plants.time], 'ascend');
    [~, orderPollinatorsIndices] = sort([pollinators.time], 'ascend');
    
    plotHelper.timePlot(plants(orderPlantIndices), pollinators(orderPollinatorsIndices), 0.75, co(3,:), co(2,:), false);
    
    xlim([0 1]);
    xlabel("Active Period");
end

function plotExampleNetwork(network)
    co = colororder();
    networkPlot = network.graphPlot(co(3,:), co(2,:), false);
    networkPlot.LineWidth = networkPlot.LineWidth * 4;
    networkPlot.MarkerSize = 10;
end

function plotRobustnessVsTheta(data, markers)

    co = colororder;
    minimumReproductionTime = unique([data.minimumReproductionTime]);
    
    plantWidth = data(1).plantWidth;
    interactionProbability = data(1).interactionProbabilities;
    numberOfPlants = data(1).numberOfPlants;
    
    for i = 1:length(data)
        data(i).meanRobustness = mean([data(i).robustnessCalculations], 'all');
    end
    
    for i = 1:length(minimumReproductionTime)
        dataSubset = data([data.minimumReproductionTime] == minimumReproductionTime(i) & [data.pollinatorWidth] > 0); %remove w=0 as these are always r=0
        scatter([dataSubset.pollinatorWidth] ./ [dataSubset.plantWidth], [dataSubset.meanRobustness], LineWidth=2, Marker=markers(i));
        hold on;
    end
    
    pollinatorWidths = unique([data.pollinatorWidth]);
    approximateRobustnesses = zeros(length(pollinatorWidths), 1);
    
    p = (plantWidth + pollinatorWidths) * interactionProbability * numberOfPlants / (1 - plantWidth);
    aggregatedRobustnessApprox = 1 - 1./p.*(1-exp(-p));
        
    for i = 1:length(approximateRobustnesses)
        approximateRobustnesses(i) = calculateApproximateRobustness(data(1).plantWidth, pollinatorWidths(i), data(1).numberOfPlants, data(1).interactionProbabilities / (1 - data(1).plantWidth));
    end

    effectivePlantDensity = data(1).numberOfPlants * data(1).interactionProbabilities * plantWidth / (1 - plantWidth);

    for i = 2:length(minimumReproductionTime)-1
        x = 0:0.01:1 / minimumReproductionTime(i);
        if minimumReproductionTime(i) < 0.5
            y = 1 - 1./(effectivePlantDensity * (1 + (1 - 2 * minimumReproductionTime(i)) .* x));
        else
            y = 1 - 1 / effectivePlantDensity * (1 + (2 * minimumReproductionTime(i) - 1) .* x);
        end
        plot(x, y, LineWidth=2, LineStyle="-", Color=co(i,:));
        hold on;
    end

    plot(pollinatorWidths / plantWidth, aggregatedRobustnessApprox, LineWidth=2, Color=co(1,:));
    hold on;
    plot(pollinatorWidths / plantWidth, approximateRobustnesses, LineWidth=2, Color=co(length(minimumReproductionTime),:));

    legendString = "a = " + minimumReproductionTime;
    
    legend(legendString, Location='southwest');
    
    
    xlabel("Interval Duration Ratio \theta");
    ylabel("Robustness R")

    ylim([0 1]);
    xlim([0 max(pollinatorWidths) / data(1).plantWidth]);
end

function plotRobustnessVsMrt(data, markers)  
    plantWidth = data(1).plantWidth;
    pollinatorWidths = unique([data.pollinatorWidth]);
        
    for i = 1:length(data)   
        data(i).robustnessMean = mean([data(i).robustnessCalculations], 'all');
    end
     
    colourOrder = get(gca,'ColorOrder');
    
    for i=1:length(pollinatorWidths)
        dataSubset = data([data.pollinatorWidth] == pollinatorWidths(i)); 
        scatter([dataSubset.minimumReproductionTime], [dataSubset.robustnessMean], LineWidth=2, Marker=markers(i), MarkerEdgeColor=colourOrder(i,:));
        hold on;   
    end
    
    effectivePlantDensity = data(1).numberOfPlants * data(1).interactionProbabilities * plantWidth / (1 - plantWidth);
    
    for i=1:length(pollinatorWidths)
        predictionMaxA = plantWidth/pollinatorWidths(i);
        a = 0:0.01:1;
        predictedR = zeros(length(a), 1);
    
        for j = 1:length(a)
            if a(j) > predictionMaxA
                break;
            end
            if a(j) < 0.5
                predictedR(j) = 1 - 1/(effectivePlantDensity * (1 + (1 - 2*a(j))*pollinatorWidths(i)/data(1).plantWidth));
            else 
                predictedR(j) = 1 - ((1 + (2*a(j) - 1)*pollinatorWidths(i)/data(1).plantWidth)) / effectivePlantDensity;
            end
        end
        plot(a(predictedR ~= 0), predictedR(predictedR ~= 0), LineWidth=2, LineStyle='-', Color=colourOrder(i,:));
        hold on;
    end
    legendString = "\theta = " + pollinatorWidths / plantWidth;
    legend(legendString, Location='southwest');
    
    xlabel("Breeding Fraction a");
    
    xlim([0 1]);
    ylim([0 1]);
end

function data = getRobustnessData( ...
    numberOfPlants, ...
    numberOfPollinators, ...
    plantWidth, ...
    interactionProbability, ...
    pollinatorWidths, ...
    minimumReproductionTimes, ...
    iterations, ...
    extinctionSequences, ...
    numberOfWorkers)

    data = struct( ...
    'interactionProbabilities', interactionProbability, ...
    'plantWidth', plantWidth, ...
    'pollinatorWidth', num2cell(repmat(pollinatorWidths, 1, length(minimumReproductionTimes))), ...
    'minimumReproductionTime', num2cell(repelem(minimumReproductionTimes, length(pollinatorWidths))), ...
    'iterations', iterations, ...
    'numberOfPlants', numberOfPlants, ...
    'numberOfPollinators', numberOfPollinators);

    plantWidthDistribution = makedist('Normal', mu = plantWidth, sigma = 0);

    dataPoints = length(data);

    parfor (i = 1:length(data), numberOfWorkers)
        robustnessCalculations = zeros(iterations, 1);
    
        pollinatorWidthDistribution = makedist('Normal', mu = data(i).pollinatorWidth, sigma = 0);
    
        networks = TemporalNetwork.createFromDistributions( ...
            plantWidthDistribution, ...
            pollinatorWidthDistribution, ...
            numberOfPlants, ...
            numberOfPollinators, ...
            true, ...
            [iterations, 1], ...
            interactionProbability, ...
            data(i).minimumReproductionTime);

        for j = 1:length(networks)
            robustnessCalculations(j) = mean(networks(j).simulateExtinctions(extinctionSequences).robustnesses);
        end        
        data(i).robustnessCalculations = robustnessCalculations;
        disp("completed: " + i + "/" + dataPoints);
    end
end

function robustness = calculateApproximateRobustness(plantWidth, pollinatorWidth, numberOfPlants, interactionProbability, breedingFraction)
arguments
    plantWidth 
    pollinatorWidth 
    numberOfPlants 
    interactionProbability 
    breedingFraction (1,1) double = 1
end
    robustness = ones(length(plantWidth), 1);

    rho = plantWidth .* interactionProbability .* numberOfPlants;
    theta = pollinatorWidth ./ plantWidth;

    if breedingFraction < 1
        if breedingFraction * theta > 1
            robustness = 0;
            warning("a\theta > 1 - no approximation given");
            return;
        end
        if breedingFraction < 0.5
            robustness = 1 - 1/(numberOfPlants * interactionProbability * plantWidth * (1 + (1 - 2 * breedingFraction) * pollinatorWidth / plantWidth));
        else 
            robustness = 1 - ((1 + (2 * breedingFraction - 1) * pollinatorWidth / plantWidth)) / (numberOfPlants * interactionProbability * plantWidth);
        end
        return;
    end

    for i=1:length(robustness)
        summation = 0;
        for j=0:floor(pollinatorWidth(i)/plantWidth(i))
            innerSum = 0;
                for l=0:j+1
                    innerSum = innerSum + ((j+1) * rho(i))^(l)/factorial(l);
                end
            %summation = summation + (j - pollinatorWidth(i)/plantWidth(i))^j/(j+1)^(j+2)*(1 - exp(- (j+1) * rho(i)) * (innerSum + (theta(i) - j) * ((j+1) * rho(i))^(j+1)/ (factorial(j+1) * (1+theta(i)))));
            summation = summation + (j - pollinatorWidth(i)/plantWidth(i))^j/(j+1)^(j+2);
        end
        robustness(i) = robustness(i) - summation * (pollinatorWidth(i) + plantWidth(i))/((numberOfPlants(i)+1) * plantWidth(i)^2 * interactionProbability(i));
    end
    % robustness = 1 - (plantWidth + pollinatorWidth) ./ ((numberOfPlants + 1) .* plantWidth.^2); 
    % 
    % isPollinatorWider = pollinatorWidth > plantWidth;
    % 
    % robustness = robustness + (pollinatorWidth.^2 - plantWidth.^2) ./ (8 .* (numberOfPlants + 1) .* plantWidth.^3) .* isPollinatorWider;
end
