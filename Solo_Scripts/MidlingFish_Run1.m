%% This code takes the mask coordinates provided by the user, creates a midline, stores both the midline coordinates in an output folder AND the overlaid midline images
%%
v = VideoReader('[insert filepath to the video you_re analysing to get its frame rate]');
framefreq=3; % How often you sampled images when creating the masks

% Get the frame rate (frames per second)
frameRate = v.FrameRate/framefreq;

% Calculate the delay time from the frame rate
delayTime =( 1 / frameRate);
%% Repository of fish files

%Your Fish Project
% Define the directory containing the text files with x and y points
txtFilesDir = '[insert filepath to your labels]';
% Define the directory containing the original images
imagesDir = '[insert filepath to your sampled images]';
% Output video file setup
outputFolder='[create an output folder first, then send files there]';
outputVideoFile = '[filepath to your output video].mp4';



%%

% Get a list of all text files in the specified directory
txtFiles = dir(fullfile(txtFilesDir, '*.txt'));
imageFiles=dir(fullfile(imagesDir,'*.jpg'));
fileIdxi=0;

% Define the directory containing the original images
imageDirectory= imagesDir;

% List all image files in the directory
imageFiles = dir(fullfile(imageDirectory, '*.jpg')); % Change the file extension if needed

% Create VideoWriter object
v = VideoWriter(outputVideoFile, 'MPEG-4');
v.FrameRate = frameRate; % Set to desired frame rate
open(v);

% Loop through each text file and corresponding image
for fileIdx = 1:numel(txtFiles)
    disp(fileIdx);
    % Read text file for coordinates
    txtFileName = txtFiles(fileIdx).name;
    txtFilePath = fullfile(txtFilesDir, txtFileName);
    fid = fopen(txtFilePath, 'r');
    data = textscan(fid, '%f');
    fclose(fid);
    f2 = data{1};

    % Separate x and y coordinates
    f = f2(2:end);
    x_d = f(1:2:end);
    y_d = f(2:2:end);

    % Load corresponding image
    imageFileName = imageFiles(fileIdx).name;
    imageFilePath = fullfile(imagesDir, imageFileName);
    originalImage = imread(imageFilePath);

    % Ensure the original image is in RGB format
    if size(originalImage, 3) == 1
        originalImage = repmat(originalImage, [1, 1, 3]);
    end

    [height, width, ~] = size(originalImage);

    % Create binary mask
    binaryMask = false(height, width);
    pixelX = round(x_d * width);
    pixelY = round(y_d * height);

    % Get the lengths of pixelX and pixelY
    lenX = length(pixelX);
    lenY = length(pixelY);
    
    % Adjust the size of pixelX and pixelY to match the smaller length
    if lenX > lenY
        pixelX(lenY+1:end) = [];  % Truncate pixelX
    elseif lenY > lenX
        pixelY(lenX+1:end) = [];  % Truncate pixelY
    end
    validIdx = pixelX > 0 & pixelX <= width & pixelY > 0 & pixelY <= height;
    binaryMask(sub2ind([height, width], pixelY(validIdx), pixelX(validIdx))) = true;
    
    % Apply the mask to the image
    pen = strel('disk', 20);  % Structuring element for morphological operations
    binaryMask = imerode(binaryMask, pen);  % Erode the mask
    binaryMask = bwpropfilt(binaryMask, 'Area', 1, 'largest');  % Keep the largest component
    binaryMask = imdilate(binaryMask, pen);  % Dilate to restore shape
    
    % Convert the binary mask into a filled mask
    binary_mask1 = poly2mask(pixelX, pixelY, height, width);
    binary_mask2 = imfill(binary_mask1, 'holes');

    % the standard skeletonization:
    imshow(bwmorph(binary_mask2,'skel',inf));
        
    % Better Skeletonization Method [Credit to Nicholas Howe: Matlab file exchange 'Better Skeletonization' function. See ReadMe]:
    imshow(bwmorph(binary_mask2>100,'skel',Inf));   %was 35 when i started
        
    
    [skr,rad] = skeleton(binary_mask2);   
    imagesc(skr);
    colormap jet
    axis image off    
    imagesc(rad)
    colormap jet
    axis image off

    %Adjust the value of 'skr>[val]' if desired
    skel = bwmorph(skr > 50,'skel',inf);
    imshow(skel)

    % anaskel returns the locations of endpoints and junction points
    [dmap,exy,jxy] = anaskel(skel);
        hold on
    plot(exy(1,:),exy(2,:),'go')
    plot(jxy(1,:),jxy(2,:),'ro')
    %End Better Skeletonization Method

    % code to sort the midline co-ordinates
    [row,col] = find(skel);              %row and col are midlines 
    A = [col row];                       %reshape the matrix. row 1 in the image is higher than row 200, so change row 200 to -200 to make sense in a plot
    points = length(col);   

    k = 2;
    distances = [];
    ordered = [A(1,:)];

    x_orig = A(1,1); y_orig = A(1,2);   %store initial point
    x = x_orig; y = y_orig;             %initialize for algorithm
    num_search = length(A(:,1));
    while length(A) > 0
        k1 = k;                         %store index of initial point

        while abs(x - A(k,1)) <= 2e5      %search until points are far from initial point in x
            dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 );     %compute distance from initial point
            distances = [distances; k, dist];                   %array of distances and point index
            k = k+1;                                            %search the next point
            if k > num_search                                   %break if at the end of the array
                break;
            end
        end

        k = k1-1;                       %search from initial point in other direction
        if k > 0                        %if not at the start of the array
            while abs(x - A(k,1)) <= 2e5 && k > 0                 %search until points are far from initial point in x && not at the start of array
                dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 ); %compute distance from initial point
                distances = [distances; k, dist];               %array of distances and point index
                k = k-1;                                        %search the next point
                if k <= 0                                       %break if at the end of the array
                    break;
                end
            end
        end

        [min_dist,index] = min(distances(:,2));                 %out of all the nearby points, which is closest?
        if min_dist < 10                                        %if not at the end of the curve
            index = distances(index,1);                         
            ordered = [ordered; A(index,:)];                    %add the nearest point to the end of the ordered array
            x = A(index,1);                                     %save this point as the new initial point for next search
            y = A(index,2);
            A = [A(1:index-1,:);A(index+1:num_search,:)];       %remove the point from the array of to-be-sorted points
            k = index;                                          %
            distances = [];                                     %clear the array of distances computed
            num_search = length(A(:,1));    %find max number of points to search
        else
            break;
        end

        if k > num_search           %if the end of the array is reached
            break;                  %end the loop
        end
    end 

    % check ordered(1) vs A(1) and A(n) and ordered(n) vs A(1) and A(n)
    if length(A) > 0    % if there are still more points to order

        endpoints = [ordered(1,:); ...
            ordered(length(ordered),:); ...
            A(1,:);...
            A(length(A(:,1)),:)];

        check_endpoints = ...       %[1,3;1,4;2,3;2,4]
        [sqrt( (endpoints(1,1) - endpoints(3,1))^2 + (endpoints(1,2) - endpoints(3,2))^2 );...  find distance between the first element of the ordered vector and first element of remaining vector
        sqrt( (endpoints(1,1) - endpoints(4,1))^2 + (endpoints(1,2) - endpoints(4,2))^2 );...   
        sqrt( (endpoints(2,1) - endpoints(3,1))^2 + (endpoints(2,2) - endpoints(3,2))^2 );...
        sqrt( (endpoints(2,1) - endpoints(4,1))^2 + (endpoints(2,2) - endpoints(4,2))^2 )];

        [min_check_endpoints,idx] = min(check_endpoints);
        switch idx
            case 1                                          %match first element of sorted vector and first element of remaining vector
                k = 1;                                      %start searching first element of remaining points
                x = endpoints(1,1); y = endpoints(1,2);     %compare distance to first element of sorted vector
            case 2                                          %match first element of sorted vector and last element of remaining vector
                k = length(A(:,1));                         %start searching last element of remaining points
                x = endpoints(1,1); y = endpoints(1,2);     %compare distance to first element of sorted vector
            case 3                                          %match last element of sorted vector and first element of remaining vector
                k = 1;                                      %start searching first element of remaining points
                x = endpoints(2,1); y = endpoints(2,2);     %compare distance to last element of sorted vector
            case 4                                          %match last element of sorted vector and last element of remaining vector
                k = length(A(:,1));                         %start searching last element of remaining points
                x = endpoints(2,1); y = endpoints(2,2);     %compare distance to last element of sorted vector
        end
        num_search = length(A(:,1));                        %find max number of points to search

        while length(A) > 0                                 %until all points have been sorted
            k1 = k;                                         %store index of initial point

            while abs(x - A(k,1)) <= 2e5                      %search until points are far from initial point in x
                dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 );     %compute distance from initial point
                distances = [distances; k, dist];                  %array of distances and point index
                k = k+1;                                            %search the next point
                if k > num_search                                   %break if at the end of the array
                    break;
                end
            end

            k = k1-1;                                               %search from initial point in other direction
            if k > 0                                                %if not at the start of the array
                while abs(x - A(k,1)) <= 2e5 && k > 0                 %search until points are far from initial point in x && not at the start of array
                    dist = sqrt( (x - A(k,1))^2 + (y - A(k,2))^2 ); %compute distance from initial point
                    distances = [distances; k, dist];               %array of distances and point index
                    k = k-1;                                        %search the next point
                    if k <= 0                                       %break if at the end of the array
                        break;  
                    end
                end
            end

            [min_dist,index] = min(distances(:,2));                 %out of all the nearby points, which is closest?
            if min_dist < 20                                        %if not at the end of the curve
                index = distances(index,1);
                if idx <= 2
                    ordered = [A(index,:);ordered];                 %add the nearest point to the start of the ordered array
                elseif idx > 2
                    ordered = [ordered; A(index,:)];                %add the nearest point to the end of the ordered array
                end
                k = index;
                x = A(index,1);                                     %save this point as the new initial point for next search
                y = A(index,2);
                A = [A(1:index-1,:);A(index+1:num_search,:)];       %remove the point from the array of to-be-sorted points

                distances = [];                                     %clear the array of distances computed
                num_search = length(A(:,1));    %find max number of points to search
                if k > num_search               %make sure not to go out of array bounds
                    k = k-1;
                end
            else
                break;
            end
        end 
    end

        if ordered(1,1) > ordered(length(ordered(:,1)),1)       %first element of ordered vector
            ordered = flip(ordered);                            %is forced to be the most negative x-coord

        elseif ordered(1,1) == ordered(length(ordered(:,1)),1)  %if endpoints are the same x-coord
            if ordered(1,2) > ordered(length(ordered(:,2)),2)   %first element of ordered vector
                ordered = flip(ordered);                        %is forced to be the most negative y-coord
            end
        end


    % Define directories to save midline coordinates
    coordsOutputFolder = fullfile(outputFolder, 'midline_coordinates');
    if ~exist(coordsOutputFolder, 'dir')
        mkdir(coordsOutputFolder);
    end
    coordsFileName = fullfile(coordsOutputFolder, [txtFileName]);
    fidCoords = fopen(coordsFileName, 'w');

    if fidCoords == -1
        error('Failed to open the file for writing: %s', coordsFileName);
    end
   

    for i = 1:length(ordered(:,1))
        fprintf(fidCoords, '%f %f\n', ordered(i,1)/width, ordered(i,2)/height); %ognl
        %fprintf(fidCoords, '%f %f\n', X_mid(i), Y_mid(i));
    end
    fclose(fidCoords);

    %Output image directory
    overlayFolder = fullfile(outputFolder, 'OverlaidMidlineImages');
    if ~exist(overlayFolder, 'dir')
        mkdir(overlayFolder);
    end

    % Dilate the skeleton for better visibility
    se = strel('disk', 1); % Adjust the size for the thickness of the centerline
    dilated_skeleton = imdilate(skel, se);

    % Convert the skeleton to an RGB image with the centerline color
    centerline_color = [1, 0, 0]; % red color for the centerline
    skel_rgb = cat(3, dilated_skeleton * centerline_color(1), ...
                      dilated_skeleton * centerline_color(2), ...
                      dilated_skeleton * centerline_color(3));

    % Ensure the skeleton is a logical array
    skel_rgb = logical(skel_rgb);

    % Overlay the skeleton on the original image
    overlay2 = originalImage;
    for c = 1:3
        channel = overlay2(:,:,c);
        channel(skel_rgb(:,:,c)) = 255 * centerline_color(c);
        overlay2(:,:,c) = channel;
    end
    overlayImageName = fullfile(overlayFolder, imageFileName); % Use the same name as the input image file
    imwrite(overlay2, overlayImageName);
    % Write frame to video
    writeVideo(v, overlay2);
end

% Close the video writer
close(v);

disp('Video creation completed.');

