
% Function to calculate errors for a segment
function totalError = calculateSegmentErrors(midlinePoints, startSegment, endSegment)
    % Find the nearest points on the midline to startSegment and endSegment
    [~, indx1] = min(pdist2(midlinePoints, startSegment));
    [~, indx2] = min(pdist2(midlinePoints, endSegment));

    % Extract the midline points between the nearest points
    midlineSegment = midlinePoints(min(indx1, indx2):max(indx1, indx2), :);

    % Calculate the vector representing the line segment
    v1 = endSegment - startSegment;

    % Calculate the coefficients of the line equation Ax + By + C = 0
    A = -v1(2);
    B = v1(1);
    C = -(A * startSegment(1) + B * startSegment(2));

    % Calculate the distance of each midline point to the line
    distances = abs(A * midlineSegment(:, 1) + B * midlineSegment(:, 2) + C) / norm([A, B]);

    % Sum all distances
    totalError = sum(distances);
end