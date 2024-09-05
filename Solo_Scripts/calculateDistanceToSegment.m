% Function to calculate the distance from a point to a line segment
function distance = calculateDistanceToSegment(point, segmentPoints)
    % Calculate the vectors representing the line segment
    v1 = segmentPoints(2, :) - segmentPoints(1, :);
    v2 = point - segmentPoints(1, :);

    % Calculate the projection of v2 onto v1
    t = dot(v2, v1) / dot(v1, v1);

    % Check if the projection is beyond the segment bounds
    if t < 0
        distance = norm(point - segmentPoints(1, :));
    elseif t > 1
        distance = norm(point - segmentPoints(2, :));
    else
        % Calculate the distance to the line segment
        projection = segmentPoints(1, :) + t * v1;
        distance = norm(point - projection);
    end
end