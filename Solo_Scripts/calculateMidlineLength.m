function totalLength = calculateMidlineLength(midlinePoints)
    totalLength = 0;
    
    for j = 2:size(midlinePoints, 1)
        x1 = midlinePoints(j-1, 1);
        y1 = midlinePoints(j-1, 2);
        x2 = midlinePoints(j, 1);
        y2 = midlinePoints(j, 2);

        % Calculate the distance between consecutive points
        distance = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        totalLength = totalLength + distance;
    end
end