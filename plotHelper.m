classdef plotHelper
    %PLOTHELPER helper class for methods useful for plotting
    methods(Static)
        function linearFitCoefficients = linearFit(x, y, includeIntercept)
            %LINEARFITCOEFFICIENTS gets coefficients for a linear regression
            %   Returns gradient and intercept for a line of best fit
            %   calulated via linear regression. If includeIntercept is 
            %   false, the line of best fit will pass through the origin
            arguments
                x {mustBeNumeric, mustBeVector}
                y {mustBeNumeric, mustBeVector}
                includeIntercept (1, 1) logical = true
            end
            if includeIntercept
                coefficients = [ones(length(x), 1), reshape(x, length(x), 1)] \ reshape(y, length(y), 1);

                linearFitCoefficients.intercept = coefficients(1);
                linearFitCoefficients.gradient = coefficients(2);
            else
                coefficients = reshape(x, length(x), 1) \ reshape(y, length(y), 1);

                linearFitCoefficients.intercept = 0;
                linearFitCoefficients.gradient = coefficients(1);
            end
        end

        function f = fillError(x, yMin, yMax, colour, alpha)
            %FILLERROR plots the lines yMax against x and yMin against x,
            %and fills the region bounded between the two lines
            arguments
                x {mustBeNumeric, mustBeVector}
                yMin {mustBeNumeric, mustBeVector}
                yMax {mustBeNumeric, mustBeVector}
                colour
                alpha (1, 1) {mustBeInRange(alpha, 0, 1, "inclusive")} = 0.5
            end

            x = reshape(x, 1, length(x));
            yMin = reshape(yMin, 1, length(yMin));
            yMax = reshape(yMax, 1, length(yMax));

            x2 = [x, flip(x)];
            inBetween = [yMin, flip(yMax)];
            f = fill(x2, inBetween, colour);
            f.FaceAlpha = alpha;
        end

        function p = timePlot(plants, pollinators, barWidth, plantColour, pollinatorColour, name)
            %TIMEPLOT plots timeline (Gantt chart) for a given set of
            %plants and pollinators
            arguments
                plants Species {mustBeVector}
                pollinators Pollinator {mustBeVector}
                barWidth (1,1) {mustBeNumeric} = 1
                plantColour = [1 0 0]
                pollinatorColour = [0 1 0]
                name (1, 1) logical = true
            end
            
            y = repelem(repmat((1:(length(plants) + length(pollinators))), 2, 1) + [-1; 1] * barWidth / 2, 2, 1);

            x = [[pollinators.startTime; pollinators.endTime; pollinators.endTime; pollinators.startTime], ...
                [plants.startTime; plants.endTime; plants.endTime; plants.startTime]];

            c = [ones(length(pollinators), 1); zeros(length(plants), 1)];

            p = patch(x,y,c);
            grid on;
            colormap([plantColour; pollinatorColour]);


            ax = gca();

            if name
                ax.YTick = 1:(length(plants) + length(pollinators));
                ax.YTickLabel = [pollinators.name plants.name];
            else
                ax.YTick = [0.5 + length(pollinators) / 2, 0.5 + length(plants) / 2 + length(pollinators)];
                ax.YTickLabel = ["Pollinators" "Plants"];
            end


            ylim([1 - barWidth, length(plants) + length(pollinators) + barWidth]);
        end
    end
end