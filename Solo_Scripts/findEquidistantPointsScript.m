 function [equidistantPoints,storeholder, before, after, totalLength] = findEquidistantPointsTrial(midlinePoints, numPoints, headsize_factor, prev_head, fileIdx)
%
    % Calculate the total length of the midline
    totalLength = calculateMidlineLength(midlinePoints);

    %Head is kept a constant size, so we must subtract its length from the length of the midline
    headlessLength=totalLength*(1-headsize_factor);
    
    %Comment/Uncomment if frame 1's midline is the right way round
    %midlinePoints=flipud(midlinePoints);      

    %hold where the head point is
    if fileIdx>1
        top_to_phead= sqrt((midlinePoints(1,1) - prev_head(1,1))^2 + (midlinePoints(1,2) - prev_head(1,2))^2);
        bottom_to_phead= sqrt((midlinePoints(end,1) - prev_head(1,1))^2 + (midlinePoints(end,2) - prev_head(1,2))^2);
        if top_to_phead>bottom_to_phead
           midlinePoints=flipud(midlinePoints); 
        end
        %add a second check, first one is imperfect:
        top_x_to_p_head= sqrt((midlinePoints(1,1) - prev_head(1,1))^2);
        top_y_to_p_head=sqrt((midlinePoints(1,2) - prev_head(1,2))^2);
        bottom_x_to_p_head= sqrt((midlinePoints(end,1) - prev_head(1,1))^2);
        bottom_y_to_p_head= sqrt((midlinePoints(end,2) - prev_head(1,2))^2);
        if top_x_to_p_head>bottom_x_to_p_head && top_y_to_p_head>bottom_y_to_p_head
            midlinePoints=flipud(midlinePoints);
        end
    end


    % Calculate the spacing between equidistant points

    spacing=headlessLength/(numPoints-1);

     % Initialize variables
    equidistantPoints = zeros(numPoints, 2);
    equidistantPoints(1, :) = midlinePoints(1, :); % Head point

    % Find equidistant points along the midline
    storeholder=zeros(numPoints,1);
    j = 2;
    new_point = 0;


    %Discard the coordinates between headlessLength and the totalLength. Change where spacing begins. Start segmenting it from spacing.
    % Skip the initial head length
    accumulatedLength = 0;
    while accumulatedLength < totalLength * headsize_factor
        x1 = midlinePoints(j-1, 1);
        y1 = midlinePoints(j-1, 2);
        x2 = midlinePoints(j, 1);
        y2 = midlinePoints(j, 2);

        dist_midline_pts = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        accumulatedLength = accumulatedLength + dist_midline_pts;
        j = j + 1;
    end

    %Hold which coordinate came before and after a given point
    before= zeros(numPoints,2);
    after= zeros(numPoints,2);

    for i = 2:numPoints
      current_length = 0;
      prev_length=0;

      while j <= size(midlinePoints,1)
		if new_point ==1
			x1 = x_new;
			y1 = y_new;
			new_point = 0;
		else
			x1 = midlinePoints(j-1, 1);
                	y1 = midlinePoints(j-1, 2);
                end
		x2 = midlinePoints(j, 1);
                y2 = midlinePoints(j, 2);

		dist_midline_pts = sqrt((x2 - x1)^2 + (y2 - y1)^2);

		prev_length = current_length;
		current_length = current_length + dist_midline_pts;

		if current_length >= spacing
			break;
		end
		j=j+1;
	end

	if current_length > spacing
	% we need to go back to prev_length (x1,y1) and add the leftover amount, not the overshoot amount!

		remaining_length = spacing - prev_length;
		x_new = x1 + (remaining_length/dist_midline_pts) * (x2-x1);
		y_new = y1 + (remaining_length/dist_midline_pts) * (y2-y1);
        after(i,1)=x2;
        after(i,2)=y2;
        before(i,1)=x1;
        before(i,2)=y1;
	else
		x_new = x2;
		y_new = y2; 
	end

	equidistantPoints(i,:) = [x_new, y_new];
	new_point = 1;
	dist = prev_length +  sqrt((x_new - x1)^2 + (y_new - y1)^2);
	storeholder(i-1) = dist;
    end
    storeholder(numPoints)=spacing;


 end


