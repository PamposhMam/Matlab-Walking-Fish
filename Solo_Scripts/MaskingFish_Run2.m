%% This code takes the masked images created in section 1 and overlays the mask YOU created over them, thus giving you an image of the fish, the mask and the midline
%%
% Output video file setup
outputVideoFile = '[output file location].mp4';
MidlineImages= '[directory to the midlined images]';

%%
% Get lists of files
txtFiles = dir(fullfile(txtFilesDir, '*.txt'));
imageFiles = dir(fullfile(MidlineImages, '*.jpg'));

% Create VideoWriter object
v = VideoWriter(outputVideoFile, 'MPEG-4');
v.FrameRate = frameRate; % Set to desired frame rate
open(v);

% Create output folder for midline coordinates and overlaid images
outputFolder = fullfile(outputFolder, 'Fish1Midline');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
overlayFolder = fullfile(outputFolder, 'OverlaidImages');
if ~exist(overlayFolder, 'dir')
    mkdir(overlayFolder);
end

% Loop through each text file and corresponding image
for fileIdx = 1:numel(txtFiles)
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

    % Save the midline coordinates to a text file
    midlineFileName = fullfile(outputFolder, txtFileName); % Use the same name as the input text file
    fidMidline = fopen(midlineFileName, 'w');
    fprintf(fidMidline, '%f\n', [x_d; y_d]);
    fclose(fidMidline);

    % Load corresponding image
    imageFileName = imageFiles(fileIdx).name;
    imageFilePath = fullfile(MidlineImages, imageFileName);
    originalImage = imread(imageFilePath);
    [height, width, ~] = size(originalImage);

    % Create binary mask
    binaryMask = false(height, width);
    pixelX = round(x_d * width);
    pixelY = round(y_d * height);
    validIdx = pixelX > 0 & pixelX <= width & pixelY > 0 & pixelY <= height;
    binaryMask(sub2ind([height, width], pixelY(validIdx), pixelX(validIdx))) = true;
    binary_mask1 = poly2mask(pixelX, pixelY, height, width);
    binary_mask2 = imfill(binary_mask1, 'holes');

    % Create skeleton
    skel = bwmorph(binary_mask2, 'skel', Inf);
    while true
        e = bwmorph(skel, 'endpoints');
        if nnz(e) > 2
            skel = skel & ~e;
        else
            break;
        end
    end

    % Create the overlay image by combining the binary mask and the skeleton
    combinedMask = binary_mask2 | skel;
    
    % Convert the combined mask to an RGB image
    overlay = uint8(combinedMask) * 255; % Convert logical mask to uint8
    colorOverlay = cat(3, overlay, zeros(size(overlay)), zeros(size(overlay))); % Red color
    
    % Create a transparent overlay image
    alpha = 0.5; % Transparency factor
    overlayImage = imoverlay_alpha(originalImage, combinedMask, [1 0 0], alpha); % Overlay with transparency

    % Save the overlaid image
    overlayImageName = fullfile(overlayFolder, imageFileName); % Use the same name as the input image file
    imwrite(overlayImage, overlayImageName);

    % Write frame to video
    writeVideo(v, overlayImage);
end

% Close the video writer
close(v);

disp('Video creation completed.');

% Function to overlay mask with color and alpha blending
function overlayImage = imoverlay_alpha(image, mask, color, alpha)
    redChannel = image(:,:,1);
    greenChannel = image(:,:,2);
    blueChannel = image(:,:,3);
    
    redChannel = uint8(double(redChannel) * (1 - alpha) + double(mask) * color(1) * 255 * alpha);
    greenChannel = uint8(double(greenChannel) * (1 - alpha) + double(mask) * color(2) * 255 * alpha);
    blueChannel = uint8(double(blueChannel) * (1 - alpha) + double(mask) * color(3) * 255 * alpha);
    
    overlayImage = cat(3, redChannel, greenChannel, blueChannel);
end
