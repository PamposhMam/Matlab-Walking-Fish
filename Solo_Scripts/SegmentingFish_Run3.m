close all;

% Code to Find the video frame rate:

v = VideoReader('[original fish video]');
SegmentNum= 5;

% Get the frame rate (frames per second)
frameRate = v.FrameRate/framefreq;
headsize_factor= 0.10; %change this to whatever percentage of the midline is just the head 

% Define the directories
contours = '[filepath to original mask labels folder]';
midline = '[filepath to midline coordinates folder]';
imageFolder = '[filepath to images folder]';
outputFolder = '[filepath to outputs folder]';
outputVideoFile = '[filepath to output video].mp4';
framefreq= 3; %Change based on the frame sampling frequency you used
extension_length_head = 5;  % Number of points to add to the head
extension_length_tail = 10;  % Number of points to add to the tail
% Parameters for extension
extend_head = true;  % Set to true to enable head extension
extend_tail = true;  % Set to true to enable tail extension

%%

finFiles = dir(fullfile(finmidline, '*.txt'));
midlineFiles=dir(fullfile(midline,'*.txt'));
conFiles=dir(fullfile(contours,'*.txt'));

% Get the frame rate (frames per second)
frameRate = v.FrameRate/framefreq;
headsize_factor= 0.10; %change this to whatever percentage of the midline is just the head 

% Define the directories
midline = finmidline;

% Create VideoWriter object
v = VideoWriter(outputVideoFile, 'MPEG-4');
v.FrameRate = frameRate;
open(v);

% Get a list of all text files
midlineFiles = dir(fullfile(midline, '*.txt'));
conFiles = dir(fullfile(contours, '*.txt'));

%Needed to prevent the head from switching place
prev_head=zeros(1,2);



% Loop through each text file
for fileIdx = 1: numel(midlineFiles)
    clear xm ym data3 head_y2 head_x2 head_x head_y Head Tail sorted_x sorted_y newmidline sorted_indices2
    clear sorted_x2 sorted_y2
    
    midlfilname = midlineFiles(fileIdx).name;
    disp(['Processing: ', midlfilname])

    % Extract frame number
    frame_number = str2double(strrep(strrep(midlfilname, 'frame_', ''), '.txt', ''));
    disp(['Frame number: ', num2str(frame_number)])
    
    % Load midline and contour data
    midFilePath = fullfile(midline, midlfilname);
    conffilespath = fullfile(contours, midlfilname);
    data2 = load(conffilespath);
    f2 = data2;
    

    [height, width, ~] = size(image)
    
    maskSize = [height, width];
    binaryMask = false(maskSize);
    
    % Extract x and y coordinates for contour
    f = f2(2:end);
    x_d = f(1:2:end);
    y_d = f(2:2:end);
    
    % Scale and round coordinates
    pixelX = round(x_d * width);
    pixelY = round(y_d * height);
    
    % Ensure x_d and y_d have the same length
    minLength = min(length(x_d), length(y_d));
    x_d = x_d(1:minLength);
    y_d = y_d(1:minLength);
    
    % Scale and round coordinates
    pixelX = round(x_d * width);
    pixelY = round(y_d * height);
    
    % Ensure coordinates are within bounds
    pixelX = max(min(pixelX, maskSize(2)), 1);
    pixelY = max(min(pixelY, maskSize(1)), 1);
    
    % Create binary mask
    for pointIdx = 1:length(pixelX)
        if pixelX(pointIdx) > 0 && pixelX(pointIdx) <= maskSize(2) && pixelY(pointIdx) > 0 && pixelY(pointIdx) <= maskSize(1)
            binaryMask(pixelY(pointIdx), pixelX(pointIdx)) = true;
        end
    end
    
    % Load midline data
    fidmid = load(midFilePath);
    data3 = fidmid;
    xm = data3(:, 1);
    ym = data3(:, 2);
    midline_p = [xm, ym];
    
    % Extend the head of the midline if required
    if extend_head
        direction_vector_head = mean(diff(midline_p(1:10, :)));
        extension_points_head = repmat(midline_p(1, :), extension_length_head, 1) - (1:extension_length_head)' * direction_vector_head;
        midline_p = [extension_points_head; midline_p];
    end

    % Extend the tail of the midline if required
    if extend_tail
        direction_vector_tail = mean(diff(midline_p(end-9:end, :)));
        extension_points_tail = repmat(midline_p(end, :), extension_length_tail, 1) + (1:extension_length_tail)' * direction_vector_tail;
        midline_p = [midline_p; extension_points_tail];
    end

    midline_s = midline_p;
    head_x_real(fileIdx)=midline_s(1,1);
    head_y_real(fileIdx)=midline_s(1,2);
    
    % Define equidistant points and segments
    numEquidistantPoints = SegmentNum +1;
    EquiPlotPoints= zeros(numEquidistantPoints,2);        %Plotting points for the segments
    
    if fileIdx==1
        keeper=zeros(fileIdx,numEquidistantPoints);
        alpha=zeros(numEquidistantPoints-1, numEquidistantPoints);
        alpha_new=zeros(numEquidistantPoints-1, numEquidistantPoints);
        alpha_wrt_head=zeros(numEquidistantPoints-1, numEquidistantPoints);
        midlength=zeros(numel(midlineFiles));
        EquiPoints3D = zeros(SegmentNum, 2, numel(midlineFiles));
    end

    
    [equidistantPoints, storeholder, before, after, midlength(:,fileIdx)] = findEquidistantPointsScript(midline_s, numEquidistantPoints, headsize_factor, prev_head, fileIdx);
    keeper(fileIdx,:)=storeholder;
    
    for i = 1:numEquidistantPoints
        EquiPlotPoints(i,:)= (equidistantPoints(i, :));
        EquiPoints3D(i, 1, fileIdx) = EquiPlotPoints(i, 1); % Store x-coordinate
        EquiPoints3D(i, 2, fileIdx) = EquiPlotPoints(i, 2); % Store y-coordinate
    end    
    
    prev_head(1,:)=EquiPlotPoints(1,:);

    figure

    % Read the corresponding image
    imageFileName = strrep(midlfilname, '.txt', '.jpg');
    imageFilePath = fullfile(imageFolder, imageFileName);
    
    if exist(imageFilePath, 'file')
        image = imread(imageFilePath);

        % Get the size of the image
        [height, width, numChannels] = size(image);
        
        % Display the pixel height
        disp(['Pixel Height: ', num2str(height)]);
            
        % Display the image
        imshow(image);
        hold on;

        axis equal
        p_x=zeros(numEquidistantPoints);
        p_y=zeros(numEquidistantPoints);

        % Number of equidistant points
        numEquidistantPoints = size(EquiPlotPoints, 1);
        
        % Preallocate an array to hold the selected colors
        selectedColors = zeros(numEquidistantPoints, 3);
        
        % Randomly select bright colors for each point
        for i = 1:numEquidistantPoints
            selectedColors(i, :) = brightColors(randi(size(brightColors, 1)), :);
        end

        hold on;
        for i = 1:numEquidistantPoints

        switch mod(i,4)
            case 0 
             plot(EquiPlotPoints(i,1) * size(image, 2), EquiPlotPoints(i,2) * size(image, 1), 'o', ...
                  'MarkerSize', 2, 'LineWidth', 2, 'Color', 'cyan');
            case 1
             plot(EquiPlotPoints(i,1) * size(image, 2), EquiPlotPoints(i,2) * size(image, 1), 'o', ...
                  'MarkerSize', 2, 'LineWidth', 2, 'Color', 'yellow');
            case 2
             plot(EquiPlotPoints(i,1) * size(image, 2), EquiPlotPoints(i,2) * size(image, 1), 'o', ...
                  'MarkerSize', 2, 'LineWidth', 2, 'Color', 'magenta');
            case 3
             plot(EquiPlotPoints(i,1) * size(image, 2), EquiPlotPoints(i,2) * size(image, 1), 'o', ...
                  'MarkerSize', 2, 'LineWidth', 2, 'Color', 'green');
        end    
            % Store the point coordinates
            p_x(i) = EquiPlotPoints(i,1)*size(image,2);
            p_y(i) = EquiPlotPoints(i,2)*size(image,1);

        end
        hold off;

            pl = line(p_x, p_y, 'Color','w', 'LineWidth', 2);
                


        hold off;

        % Save the segmented image
        saveas(gcf, fullfile(outputFolder, imageFileName));
        close(gcf);
        % Write frame to video
        writeVideo(v, imread(fullfile(outputFolder, imageFileName)));
    else
        disp(['Image file not found: ', imageFilePath]);
    end

    %Angle calculation
    for i = 1:numEquidistantPoints-2
        theta2=  atan2(equidistantPoints(i+2,2) - equidistantPoints(i+1,2), equidistantPoints(i+2,1) - equidistantPoints(i+1,1));
        theta1= atan2(equidistantPoints(i+1,2) - equidistantPoints(i,2), equidistantPoints(i+1,1) - equidistantPoints(i,1));
        alpha(fileIdx, i)=theta2-theta1;
        if alpha(fileIdx,i)>pi
              alpha(fileIdx,i)=alpha(fileIdx,i)-2*pi;
              elseif alpha(fileIdx,i)<-pi
                     alpha(fileIdx,i)=alpha(fileIdx,i)+2*pi;
        end
    end


end

% Close the video file
close(v);
disp('Processing complete');

