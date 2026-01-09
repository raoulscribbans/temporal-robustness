% plots an extinction curve for a specified set of 4 plants and 2
% pollinators, for various values of a. Values are hard-coded and
% correspond to figure 1 of the main text.

extinctionSequences = 1000;

plants = initialisePlants();
pollinators = initialisePollinators();

adjacencyMatrix = [1 0; 1 0; 1 1; 1 1];

reprodutionFractions = [0 0.5 0.75 1];

lineStyles = ["-" "--" "-." ":"];

plotExtincitonCurve(plants, pollinators, reprodutionFractions, adjacencyMatrix, lineStyles, extinctionSequences)

function plants = initialisePlants()    
    plantA = Species();
    plantA.time  = 1/16;
    plantA.intervalWidth = 1/8;
  
    plantB = Species();
    plantB.time  = 0.75;
    plantB.intervalWidth = 0.5;

    plantC = Species();
    plantC.time  = 5/8;
    plantC.intervalWidth = 3/4;
    
    plantD = Species();
    plantD.time  = 3/8;
    plantD.intervalWidth = 3/4;

    plants = [plantA, plantB, plantC, plantD];
end

function pollinators = initialisePollinators()
    
    pollinator1 = Pollinator();
    pollinator1.time = 0.5;
    pollinator1.intervalWidth = 1;
    
    pollinator2 = Pollinator();
    pollinator2.time = 5/8;
    pollinator2.intervalWidth = 1/2;

    pollinators = [pollinator1, pollinator2];
end

function plotExtincitonCurve(plants, pollinators, reprodutionFractions, adjacencyMatrix, lineStyles, extinctionSequences)

    figure();
    
    for i = 1:length(reprodutionFractions)
        network = TemporalNetwork(plants, pollinators, 0, reprodutionFractions(i)).setInteractions(adjacencyMatrix);
        robustnessData = network.simulateExtinctions(extinctionSequences);
        plot(0:1/length(plants):1, robustnessData.meanPollinatorsRemaining / length(pollinators), LineWidth=3, LineStyle=lineStyles(mod(i - 1, length(lineStyles)) + 1));
        hold on;
    end
    
    hold off;
    
    fontsize(18, "points")
    
    xlabel("Plants Removed", FontSize = 22);
    ylabel("Pollinators Remaining", FontSize = 22);
    
    set(gcf,'Position',[0 0 400 400])
    xlim([0 1]);
    ylim([0 1]);
    xticks(0:1/length(plants):1)
    legend("a = " + reprodutionFractions, Location="southwest");
end