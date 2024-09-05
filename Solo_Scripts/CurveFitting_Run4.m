alphas = alpha(:, 1:end-1);
p = cell(size(alphas,2),2);
% Verify dimensions
[numPoints, numSets] = size(alphas);

% Define the sine function for fitting
ft = fittype(@(A, B, C, x) A * sin(B * x + C), ...
    'independent', 'x', 'dependent', 'y');

% Preallocate storage for fit results and goodness of fit
fitresultssine(numSets) = struct('fitresult', [], 'gof', []); % Struct to store fit results
fitresultsline(numSets) = struct('fitresult', [], 'gof', []); % Struct to store fit results
A = zeros(numSets, 1); % Preallocate for parameter A
B = zeros(numSets, 1); % Preallocate for parameter B
C = zeros(numSets, 1); % Preallocate for parameter C

% Define the number of subplots per figure
plotsPerFigure = 6;
numFigures = ceil(numSets / plotsPerFigure);

% Perform the fit for each set of alphas
for fig = 1:numFigures
    figure;
    for j = 1:plotsPerFigure
        index = (fig-1)*plotsPerFigure + j;
        if index > numSets
            break;
        end
        
        % Compute x-values as (2 * pi * (1:numPoints)' * frequencies(j)) / frameRate
        x = (1:numPoints) / (frameRate); % Compute x-values for the current 

        % Ensure x is a column vector
        x = x(:);

        % Ensure alphas(:, index) is a column vector
        y = alphas(:, index);

        % Check dimensions
        if length(x) ~= length(y)
            error('Dimensions of x and y do not match.');
        end

        % Fit the sine function to each set of y-values in alphas
        [fitresultsine, gofsine] = fit(x, y, 'sin1');
        [fitresultline, gofline] = fit(x, y, 'poly1');

        % Store results
        fitresults(index).fitresult = fitresultsine;
        fitresults(index).gof = gofsine;
        fitresultsline(index).fitresult = fitresultline;
        fitresultsline(index).gof = gofline;
        p{(fig-1)*(plotsPerFigure)+j,1} = fitresults;
        p{(fig-1)*(plotsPerFigure)+j,2} = fitresultsline;
        
        % Extract and store the parameters
        A(index) = fitresultsine.a1;
        B(index) = fitresultsine.b1;
        C(index) = fitresultsine.c1;

        % Create subplot
        subplot(2, 3, j); % 2 rows, 3 columns
        plot(x, y, 'bo'); % Original data
        hold on;
        plot(x, feval(fitresultsine, x), 'r-'); % Fitted sine curve
        xlabel('time/s'); % Label the x-axis appropriately
        ylabel('Angle/rad');
        legend('Data', 'Fitted Sine Curve');
        title(['Sine Curve Fit for Set ', num2str(index)]);
    end
end

% Display or use the extracted parameters
disp('Extracted Parameters:');
disp(table((1:numSets)', A, B, C, 'VariableNames', {'Set', 'A', 'B', 'C'}));
%% Plot the r-stats of every point


% Initialize vectors to store indices
lessthan = zeros(numSets-1, 1);
greaterthan = zeros(numSets-1, 1);

% Initialize vector to store R-squared values
r_vals = zeros(numSets-1, 1);

% Initialize vector to store R-squared values for plotting
rsquaredframes = zeros(numSets-1, 1);

% Loop over each set (excluding the last one)
for i = 1:numSets-1
    % Extract the relevant fit results
    temp = p{numSets-1, 1};
    
    % Store the R-squared value
    rsquaredframes(i) = temp(1, i).gof.rsquare;
    
    % Classify based on the R-squared value
    if rsquaredframes(i) < 0.1
        lessthan(i, 1) = i;
    else 
        greaterthan(i, 1) = i;
    end   
    
    % Display the current index and R-squared value
    disp(i);
    disp(temp(1, i).gof.rsquare);
    
    % Store the R-squared value for plotting
    r_vals(i) = temp(1, i).gof.rsquare;
end

% Prepare to plot all points at once
figure;
scatter(1:numSets-1, r_vals, 'filled','r'); % Plot all R-squared values
xlabel('Angle Number', FontSize=15); % Label the x-axis appropriately
ylabel('R2-value', FontSize= 15);
title(['Sine Curve Fit R2-stats ', num2str(length(alpha(1,:))-1)]);
grid on
grid minor
disp('Greater than 0.1:');
disp(greaterthan(greaterthan > 0)); % Only display non-zero indices
disp('Less than 0.1:');
disp(lessthan(lessthan > 0)); % Only display non-zero indices

