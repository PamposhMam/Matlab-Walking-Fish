function [Distances] = ErrorDiagnostics(midlinePoints, equidistantPoints, numPoints)

% Segment coordinates
x_eq = equidistantPoints(:,1);
y_eq = equidistantPoints(:,2);

% Midline Points
x_mid = midlinePoints(:,1);
y_mid = midlinePoints(:,2);

% Initialize distances array
Distances = zeros(numPoints, 1);

% Calculate the total length of the midline
totalLength = calculateMidlineLength(midlinePoints);

% Calculate the spacing between equidistant points
spacing = totalLength / (numPoints - 1);

% Variables for distance calculation
dist = 0;
dist2 = 0;
j = 1;
juststarted = 1;

% Loop over each equidistant point
for i = 1:numPoints
    distthresh = i * spacing;  % Find the next Target
    
    % Cumulative midline distance
    while dist < distthresh && j < numel(x_mid)   % while we haven't exceeded our threshold distance
        if juststarted == 1  % first time, use x_eq and y_eq
            dist = dist + sqrt((x_eq(i) - x_mid(j))^2 + (y_eq(i) - y_mid(j))^2);
            juststarted = 0;
        else
            % calculate distance along midline
            dist = dist + sqrt((x_mid(j) - x_mid(j-1))^2 + (y_mid(j) - y_mid(j-1))^2);
        end
        j = j + 1;  % update j 
    end
    
    % Calculate the distance from the current equidistant point to the last midline point
    if j <= numel(x_mid)
        dist = dist + sqrt((x_eq(i) - x_mid(j-1))^2 + (y_eq(i) - y_mid(j-1))^2);
    end

    % Store the calculated distance in the Distances array
    Distances(i) = dist - dist2;  % calculating segment distance
    dist2 = dist;  % update dist2 for next iteration
    juststarted = 1;  % reset juststarted for next equidistant point
end

end

