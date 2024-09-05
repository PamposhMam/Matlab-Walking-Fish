%Script to compare the data from the simulation to the data from the original video

%Once done, save workspace with this: save('NAMEHERE!.mat', '-v7.3');

% Fish Data
fish_length_m=10.9/100; %fish length in m
fish_mass_kg=7.2/1000; %fish weight in kgs


%% Run before to initialize things if haven't already
t_diff=1/frameRate;
time_diff=t_diff;

% Define the dimensions of the image
image_height_px = height; % Height of the image in pixels
image_width_px = width; % Width of the image in pixels

%work out the length to pixel conversion ratio
fish_length_pixels=max(midlength(:,1));
m_per_px_length = fish_length_m / fish_length_pixels; % m per pixel for the length
m_per_px_width = fish_length_m / fish_length_pixels; % m per pixel for the width

% Extract head position (first point of midline_s) in each frame
head_positions_x =head_x_real'; % x-coordinates of the head
head_positions_y =head_y_real'; % y-coordinates of the head

% Convert pixel coordinates to real-world coordinates (m)
head_positions_x_m = head_positions_x * m_per_px_width;
head_positions_y_m = head_positions_y * m_per_px_length;


% Time vector for plotting
time_vector = (0:length(head_positions_x_m)-1) * time_diff;
j=0;
tim_interval=0.02;
for i=1:length(time_vector)
    if  rem( time_vector(i) , tim_interval) ==0
        j=j+1;
        interval_head_pos_x(j)=head_positions_x_m(i);
        interval_head_pos_y(j)=head_positions_y_m(i);
        time_interval_seri(j)=time_vector(i);
    end
end 


%find the delta x and the delta y, since the sim takes the head position at t=0 as the coordinate (0,0)
x_delta= zeros(fileIdx,1);
y_delta= zeros(fileIdx,1);
for i=1:fileIdx
    x_delta(i)=head_positions_x_m(i)-head_positions_x_m(1);
    y_delta(i)=head_positions_y_m(i)-head_positions_y_m(1);
end

frame_number=zeros(numel(midlineFiles));
for i=1:numel(midlineFiles)
    frame_number(i)=i;
end

[a,b,c]=size(EquiPoints3D);

% Initialize an array to store the absolute differences
absDiffs = zeros(a,1,c);
x_Diffs=zeros(size(absDiffs));
y_Diffs=zeros(size(absDiffs));


% Reference point: first x, y pair from the first slice
x_ref = EquiPoints3D(1, 1, 1); %EquiPoints3D(1, 1, i) for dynamic
y_ref = EquiPoints3D(1, 2, 1); %EquiPoints3D(1, 2, i) for dynamic

% Loop through each slice
for i = 1:size(EquiPoints3D, 3)
    % Loop through each point in the current slice
    for j = 1:size(EquiPoints3D, 1)
        % Calculate the overall distance wrt initial x,y headpoint (i.e. from frame 1)
        absDiffs(j, 1, i) = sqrt((EquiPoints3D(j, 1, i) - x_ref)^2+(EquiPoints3D(j, 2, i) - y_ref)^2);
        x_Diffs(j,1,i)=(EquiPoints3D(j,1,i)-x_ref);
        y_Diffs(j,1,i)=(EquiPoints3D(j,2,i)-y_ref);
    end
end





%%
%Normalization: Measure all of these variables wrt their own bodylength
fishwhole=(373.5)/1000;     %the divide by 1k is bc the transform sensor outputs in metres and these measurements were in mm

x_position_sim= out.out_x1.signals.values;
y_position_sim= out.out_y1.signals.values;

x_position_sim=x_position_sim/fishwhole;
y_position_sim=y_position_sim/fishwhole;

x_delta=x_delta/fish_length_m;
y_delta=y_delta/fish_length_m;
head_positions_x_m=head_positions_x_m/fish_length_m;
head_positions_y_m=head_positions_y_m/fish_length_m;

%% Accounting for the rotation in the simulation data as well as the flipping of co-ordinates:
% Only uncomment if needed!


% %If flipped about the x axis
% % x_delta=-x_delta; 
% % head_positions_x_m=-head_positions_x_m;
% % x_Diffs=-x_Diffs;
 
% %If flipped about the y axis; %vid1
% y_delta=-y_delta;
% head_positions_y_m=-head_positions_y_m;
% y_Diffs=-y_Diffs;
%% Head Posn, timeless

% Downsample the second dataset
downsample_factor = length(x_position_sim)/numel(midlineFiles);
x_position_sim_ds = x_position_sim(1:downsample_factor:end);
y_position_sim_ds = y_position_sim(1:downsample_factor:end);

% Plot head position from original video data
figure;
%plot(head_positions_x_m, head_positions_y_m, 'b-', 'LineWidth', 1.5);
plot(x_delta, y_delta,'r-', 'LineWidth', 1.5);
hold on;

% Plot Simscape data (x/y positions) on the same figure
plot(x_position_sim_ds, y_position_sim_ds, 'b-', 'LineWidth', 1.5);

% Formatting the figure
axis equal;


% Add titles and labels
title([num2str(length(alpha(1,:))-1),' Segments: Head Position in Tank']);
xlabel('x\_direction/bodylengths');
ylabel('y\_direction/bodylengths');
grid on;

% Add legend to distinguish between the plots
legend('Video Data', 'Simscape Data');

hold off;

% Calculate the axis limits based on the largest bounds of either dataset
x_min = min([min(x_delta), min(x_position_sim_ds)]);
x_max = max([max(x_delta), max(x_position_sim_ds)]);
y_min = min([min(y_delta), min(y_position_sim_ds)]);
y_max = max([max(y_delta), max(y_position_sim_ds)]);

% Ensure that the axis limits are the same for both subplots
x_limits = [x_min, x_max];
y_limits = [y_min, y_max];

% Create a new figure
figure;

% First subplot: plot the first scatter plot with star markers
subplot(1, 2, 1); % Create a subplot with 1 row and 2 columns, and select the 1st subplot
scatter(x_delta, y_delta, 75, 'b', '*', 'LineWidth', 1.5); % '75' is the marker size
title('Delta Positions');
xlabel('x\_direction/bodylengths');
ylabel('y\_direction/bodylengths');
xlim(x_limits);
ylim(y_limits);
grid on;

% Second subplot: plot the downsampled scatter plot with star markers
subplot(1, 2, 2); % Select the 2nd subplot
scatter(x_position_sim_ds, y_position_sim_ds, 75, 'r', '*', 'LineWidth', 1.5); % '75' is the marker size
title('Simscape Positions (Downsampled)');
xlabel('x\_direction/bodylengths');
ylabel('y\_direction/bodylengths');
xlim(x_limits);
ylim(y_limits);
grid on;








%% Dist_Midline

% % Number of segments
numSegments = numEquidistantPoints - 1;

plotsPerFigure = 6; %Number of plots per figure
rows = 2; % Number of rows per figure
cols = 3; % Number of columns per figure



% %Create a error holding var to hold the error for each segment
dist_errors=zeros(numSegments,1);
downsample_factor = length(x_position_sim)/numel(midlineFiles);
time_vector_sim_other=out.out_x3.time;



% % Time vector for the Simscape data
time_vector_sim = out.seg1.time; % Assuming the time vector is stored as 'time' in 'out'

% % Define your custom x-axis
x_axis_custom = frame_number * t_diff; % Assuming frame_number and t_diff are defined

% Loop through each segment and plot
for segIdx = 1:numSegments
    % Determine which figure to use
    if mod(segIdx-1, plotsPerFigure) == 0
        figure;
    end

    % Extract Simscape data for the current segment, BELIEVE ME it's the dist
    dist_position = abs(out.(['seg', num2str(segIdx)]).signals.values)/fishwhole;
    dist_position_ds=dist_position(1:downsample_factor:end);

    % Extract corresponding equidistant point from EquiPlotPoints
    equi_x_position = EquiPlotPoints(segIdx, 1);  % X-coordinate from EquiPlotPoints

    % Extract corresponding absDiff for the current segment
    abs_diff_segment = absDiffs(segIdx, 1, :);
    abs_diff_segment = reshape(abs_diff_segment, [1, c]); % Reshape to a 1D vector for plotting


    %Calculate Error
    for i=1:numel(abs_diff_segment)
        dist_errors(segIdx)= dist_errors(segIdx) - (abs_diff_segment(i)-dist_position_ds(i));
    end
    % hold on;

    % Create a subplot for each segment, arranged in a 2x3 grid
    subplot(rows, cols, mod(segIdx-1, plotsPerFigure) + 1);

    % Plot Simscape data
    plot(time_vector_sim, dist_position, 'b-', 'LineWidth', 1.5); % Simscape data in blue
    hold on;

    % Plot absDiffs against the custom x-axis
    plot(x_axis_custom, abs_diff_segment, 'r--', 'LineWidth', 1.5); % absDiffs in red

    % Add titles and labels
    title(['Segment ', num2str(segIdx), ': Simulated Midline Distance vs Video Distance']);
    xlabel('Time/s'); % or another appropriate label based on what x_axis_custom represents
    ylabel('Position/body-lengths');
    %legend('Simscape Data', 'Video Data');
    grid on;
    ylim([-inf inf]);

    hold off;
    sgtitle(['Video 1, ' num2str(numSegments),' Segments']);
end



%% X_Axis Midline Discrepancy
segIdx=1;

downsample_factor = length(x_position_sim)/numel(midlineFiles);
time_vector_sim_other=out.out_x3.time;
y_position_sim_ds = y_position_sim(1:downsample_factor:end);

%Create a error holding var to hold the error for each segment
x_errors=zeros(numSegments,1);

%For just the x_axis throwoff
for segIdx = 1:numSegments
    % Determine which figure to use
    if mod(segIdx-1, plotsPerFigure) == 0
        figure;
    end

    % Extract Simscape data for the current segment, THEN DOWNSAMPLE
    sim_x_position = out.(['out_x', num2str(segIdx)]).signals.values/(fishwhole);
    sim_x_position_ds = sim_x_position(1:downsample_factor:end);
    time_vector_sim_other_ds=time_vector_sim_other(1:downsample_factor:end);
    % Extract corresponding equidistant point from EquiPlotPoints
    equi_x_position = EquiPlotPoints(segIdx, 1);  % X-coordinate from EquiPlotPoints, not sure why this is here, please check!!!

    % Extract corresponding absDiff for the current segment
    x_diff_segment = x_Diffs(segIdx, 1, :);       %relative difference from head posn at t=0
    x_diff_segment = reshape(abs_diff_segment, [1, c]); % Reshape to a 1D vector for plotting

    % Create a subplot for each segment, arranged in a 2x3 grid
    subplot(rows, cols, mod(segIdx-1, plotsPerFigure) + 1);

    %Calculate Error
    for i=1:numel(x_diff_segment)
        x_errors(segIdx)= x_errors(segIdx) + abs(x_diff_segment(i)-sim_x_position_ds(i));
    end

    % Plot Simscape data
    plot(time_vector_sim_other_ds, sim_x_position_ds, 'b-', 'LineWidth', 1.5); % Simscape data in blue
    hold on;

    % Plot absDiffs against the custom x-axis
    plot(x_axis_custom, abs_diff_segment, 'r--', 'LineWidth', 1.5); % absDiffs in red

    % Add titles and labels
    title(['Segment ', num2str(segIdx), ': Simulated Midline X\_Posn vs Video X\_Posn']);
    xlabel('Time/s'); % or another appropriate label based on what x_axis_custom represents
    ylabel('Position/body-lengths');
    legend('Simscape Data', 'Video Data');
    grid on;
    ylim([-inf inf]);
    hold off;

end


%% Y_Axis Midline Discrepancy
segIdx=1;

%Create a error holding var to hold the error for each segment
y_errors=zeros(numSegments,1);
downsample_factor = length(x_position_sim)/numel(midlineFiles);
time_vector_sim_other=out.out_x3.time;


%For just the y_axis throwoff
for segIdx = 1:numSegments
    % Determine which figure to use
    if mod(segIdx-1, plotsPerFigure) == 0
        figure;
    end

    % Extract Simscape data for the current segment
    sim_y_position = out.(['out_y', num2str(segIdx)]).signals.values/(fishwhole);
    sim_y_position_ds = sim_y_position(1:downsample_factor:end);
    time_vector_sim_other_ds=time_vector_sim_other(1:downsample_factor:end);

    % Extract corresponding absDiff for the current segment
    y_diff_segment = y_Diffs(segIdx, 1, :);
    y_diff_segment = reshape(abs_diff_segment, [1, c]); % Reshape to a 1D vector for plotting

    % Create a subplot for each segment, arranged in a 2x3 grid
    subplot(rows, cols, mod(segIdx-1, plotsPerFigure) + 1);

    % Plot Simscape data
    plot(time_vector_sim_other_ds, sim_y_position_ds, 'b-', 'LineWidth', 1.5); % Simscape data in blue

       %Calculate Error
    for i=1:numel(y_diff_segment)
        y_errors(segIdx)= y_errors(segIdx) + abs(y_diff_segment(i)-sim_y_position_ds(i));
    end
    hold on;

    % Plot absDiffs against the custom x-axis
    plot(x_axis_custom, y_diff_segment, 'r--', 'LineWidth', 1.5); % absDiffs in red


    % Add titles and labels
    ylim([-inf inf]);
    title(['Segment ', num2str(segIdx), ': Simulated Midline Y\_Posn vs Video Y\_Posn']);
    xlabel('Time/s'); % or another appropriate label based on what x_axis_custom represents
    ylabel('Position/body-lengths');
    legend('Simscape Data', 'Video Data');
    grid on;


    hold off;
end

%% Is this segment in contact with the ground plane? (Simulation Data Analysis)
j=1;
numSegments=length(alpha(1,:))-1;

segInCont=zeros(length(out.out_y1.signals.values),numSegments); %Create a variable to store which segments are in contact with the ground
time_vector_sim_other=out.out_y3.time;  %Create a time vector
time_vector_sim_other_ds=time_vector_sim_other(1:downsample_factor:end); %Downsample to get the data shown at the times we have images for

for j=1:numSegments
    sim_force = out.(['f_bot_', num2str(j)]).signals.values;
    for i=1:length(out.out_y3.time)
        if sim_force(i)~=0
        segInCont(i,j)=j;  %If j's simulation force isn't equal to 0 at time i, then the j'th segment is in contact with the ground
        end
    end
end


%Plot the normal data
figure;
for i=1:numSegments

    scatter (time_vector_sim_other, segInCont(:,i)); % Simscape data in blue
    hold on;
    % Add titles and labels
    title(['Segment in contact with the ground']);
    xlabel('Time (s)'); % or another appropriate label based on what x_axis_custom represents
    ylabel('Segment');
    grid on;
    ylim([1 numSegments]);
    xlim([-inf inf]);

end

%Then add a second plot, downsampled so it matches the number of frames in the video
j=1;
figure;
for i=1:numSegments

    current=segInCont(:,i);
    segInCont_ds=current(1:downsample_factor:end);
    scatter (time_vector_sim_other_ds, segInCont_ds); % Simscape data in blue
    hold on;
    % Add titles and labels
    title(['Segment in contact with the ground']);
    xlabel('Time (s)'); % or another appropriate label based on what x_axis_custom represents
    ylabel('Segment)');
    grid on;
    ylim([1 numSegments]);
xlim([-inf inf]);
end
%% Joint Energy (Simulation Data Analysis)

%Plot the energy of the joints

SegIdx = 1;
rows=2;
cols=3;

for segIdx = 1:NumOfSegs
    % Determine which figure to use
    if mod(segIdx-1, plotsPerFigure) == 0
        figure;
    end

    % Extract Simscape data for the current segment
    torque = out.(['out_torque', num2str(segIdx)]).signals.values;
    angdisp = A(segIdx) * sin(B(segIdx) * out.out_torque1.time + C(segIdx));
    rotnenergy = torque .* angdisp; % Element-wise multiplication, from E=T.q

    % Create a subplot for each segment, arranged in a 2x3 grid
    subplot(rows, cols, mod(segIdx-1, plotsPerFigure) + 1);

    % Plot the rotational energy over time
    plot(out.out_torque1.time, rotnenergy, 'b-', 'LineWidth', 1.5);
    hold on;

    % Add titles and labels
    title(['Segment ', num2str(segIdx), ': Rotational Energy vs Time']);
    xlabel('Time (s)');
    ylabel('Rotational Energy (Joules)');
    legend('Rotational Energy');
    grid on;
    ylim([-inf inf]);

    hold off;
end

%%
disp("Done running!");




