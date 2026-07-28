function [out] = ais_auto(cell, threshold, smooth, chosen)
% AIS_AUTO - Automated detection and measurement of Axon Initial Segment
%
% Inputs:
%   cell      - Cell/file name prefix
%   threshold - Threshold value (0 for automatic Otsu method)
%   smooth    - Smoothing method ('gaussian' or other)
%   chosen    - Not used in current implementation
%
% Output:
%   out - 0 if analysis rejected, otherwise AIS measurements
%
% Based on Grubb & Burrone 2010 methodology

close all

% ============ Constants ============
PIXEL_TO_MICRON = 0.161;  % Microns per pixel
INTENSITY_FRACTION = 0.33; % Fraction of max fluorescence for AIS start/end
MIN_COMPONENT_SIZE = 70;   % Minimum pixels for connected components
SMOOTH_KERNEL_SIZE = 20;   % Gaussian smoothing kernel size
SMOOTH_KERNEL_SD = 2;      % Gaussian smoothing standard deviation
MORPH_ELEMENT_SIZE = 3;    % Size for morphological operations
MAX_AIS_LENGTH = 300;      % Maximum expected AIS length in pixels
MAX_BRANCH_ITERATIONS = 10; % Maximum iterations for branch removal
SPLINE_SMOOTHING = 0.3;    % Spline smoothing parameter
SLIDING_WINDOW_WIDTH = 1.5; % Width in microns for sliding mean

% ============ Load and Display Image ============
disp(cell);
filename = [cell '.tif - Processed method 2.5.tif'];
img = imread(filename);

% Extract and display channel
ch1 = img(:,:,1);
ch1_processed = ch1;  % Keep a copy of the processed image
figure(1)
imagesc(ch1)
colormap(gray)
axis square
hold on

% ============ Image Preprocessing ============
% Convert to grayscale and smooth
img_gray = mat2gray(ch1);

if strcmp(smooth, 'gaussian')
    img_smooth = imfilter(img_gray, fspecial('gaussian', [SMOOTH_KERNEL_SIZE SMOOTH_KERNEL_SIZE], SMOOTH_KERNEL_SD));
else
    img_smooth = img_gray; % No smoothing if not gaussian
end

figure(2)
imagesc(img_smooth)
colormap(gray)
axis square
title('Smoothed Image')

% ============ Thresholding and Region Selection ============
if threshold == 0
    threshold = graythresh(img_smooth); % Otsu's method
end

% Interactive thresholding loop
region_selected = false;
loop_count = 0;

while ~region_selected
    % Threshold image
    img_binary = im2bw(img_smooth, threshold);

    figure(3)
    imagesc(img_binary)
    colormap(gray)
    axis square
    title('Thresholded Image')

    % Morphological operations
    se_open = strel('square', MORPH_ELEMENT_SIZE);
    se_close = strel('square', MORPH_ELEMENT_SIZE);
    img_opened = imopen(img_binary, se_open);
    img_closed = imclose(img_opened, se_close);

    figure(4)
    imagesc(img_closed)
    colormap(gray)
    axis square
    title('Morphologically Processed')

    % Find and filter connected components
    cc = bwconncomp(img_closed);
    pixel_counts = cellfun(@numel, cc.PixelIdxList);
    small_components = find(pixel_counts < MIN_COMPONENT_SIZE);
    cc.PixelIdxList(small_components) = [];
    cc.NumObjects = length(cc.PixelIdxList);

    % User selects AIS region
    figure(1)
    disp(' ')
    disp('Zoom, then click near AIS start point')
    zoom on;
    pause;
    zoom off;
    [x_click, y_click] = ginput(1);
    x_click = round(x_click);
    y_click = round(y_click);

    % Find selected component
    label_matrix = labelmatrix(cc);
    selected_label = label_matrix(y_click, x_click);

    % If clicked on background, find nearest component
    if selected_label == 0
        binary_mask = imbinarize(label_matrix, 0.0000001);
        [~, idx_nearest] = bwdist(binary_mask);
        nearest_pixel = idx_nearest(y_click, x_click);
        selected_label = label_matrix(nearest_pixel);
    end

    ais_mask = (label_matrix == selected_label);

    figure(5)
    imagesc(ais_mask)
    colormap(gray)
    axis square
    title('Selected AIS Region')

    % Check if maximum pixel is in selected region
    if loop_count == 0
        img_double = im2double(ch1);
        max_global = max(img_double(:));
        max_ais = max(max(ais_mask .* img_double));

        if max_global == max_ais
            region_selected = true;
        else
            disp('Max pixel not in AIS region. Rethresholding...')
            threshold = (max_ais / max_global) * threshold;
            loop_count = 1;
        end
    else
        region_selected = true;
    end
end

% ============ Skeletonization ============
ais_skeleton = bwmorph(ais_mask, 'thin', Inf);

figure(6)
imagesc(ais_skeleton)
colormap(gray)
axis square
title('AIS Skeleton')

% ============ Trace Skeleton ============
% Find starting point (closest to user click)
[~, idx_dist] = bwdist(ais_skeleton);
closest_pixel = idx_dist(y_click, x_click);
[start_x, start_y] = ind2sub(size(idx_dist), closest_pixel);

% Trace skeleton and handle branches
[x_trace, y_trace] = traceSkeleton(ais_skeleton, start_x, start_y, MAX_AIS_LENGTH, MAX_BRANCH_ITERATIONS);

% ============ Display Traced Path ============
figure(7)
imagesc(ch1_processed)  % Use the processed image to match what user selected
colormap(gray)
axis square
title('Traced AIS Path')
hold on
plot(y_trace, x_trace, '.r')

% ============ Fit Spline and Calculate Profile ============
% Prepare coordinates
xy = double([y_trace; x_trace]);
t = 1:length(xy);
ts = [1:0.1:length(xy), t(end)];

% Fit spline
xy_spline = csaps(t, xy, SPLINE_SMOOTHING, ts);

% Calculate distances along spline
distances = zeros(size(xy_spline, 2), 1);
for i = 2:length(xy_spline)
    distances(i) = sqrt((xy_spline(1,i) - xy_spline(1,i-1))^2 + ...
                       (xy_spline(2,i) - xy_spline(2,i-1))^2);
end
cumulative_dist = cumsum(distances);
dist_microns = cumulative_dist * PIXEL_TO_MICRON;

% Plot spline with consistent coordinates
plot(xy_spline(1,:), xy_spline(2,:), '-y', 'LineWidth', 2)

% ============ User Validation ============
[~, ~, button] = ginput(1);
button = char(button);

if strcmp(button, 'n')
    out = 0;
    return;
end

% ============ Extract Unique Pixel Coordinates ============
[x_pixels, y_pixels, dist_along_axis] = extractUniquePixels(xy_spline, cumulative_dist, PIXEL_TO_MICRON);

% ============ Save Coordinates ============
saveCoordinates(cell, x_pixels, y_pixels);

% Reload coordinates to ensure consistency (mimicking original behavior)
savexy = [cell '_xy.txt'];
fid = fopen(savexy, 'r');
n_points = fscanf(fid, '%f', 1);
coords = fscanf(fid, '%f %f', [2, n_points]);
fclose(fid);
x_pix = coords(1,:);
y_pix = coords(2,:);

% ============ Analyze Fluorescence Profile ============
disp(cell)
close all

% Load the ORIGINAL (non-processed) image for fluorescence analysis
ch1_pic = imread([cell '.tif']);

% Display both images side by side for debugging
figure(8)
subplot(1,2,1)
imagesc(ch1_processed)
colormap(gray)
title('Processed (used for detection)')
hold on
plot(x_pix, y_pix, '-r', 'LineWidth', 2)
hold off

subplot(1,2,2)
imagesc(ch1_pic)
colormap(gray)
title('Original (used for measurement)')
hold on
plot(x_pix, y_pix, '-r', 'LineWidth', 2)
hold off

% Now show the measurement image
figure(1)
imagesc(ch1_pic)
colormap(gray)
axis square
title('Ch1')
hold on

% Plot the traced path on the original image
plot(x_pix, y_pix, '-b')

% Extract fluorescence intensities from the ORIGINAL image
[intensity_raw, intensity_smooth, intensity_sliding] = extractIntensityProfile(ch1_pic, x_pix, y_pix, PIXEL_TO_MICRON);

% Normalize intensities
norm_intensity = (intensity_sliding - min(intensity_sliding)) / (max(intensity_sliding) - min(intensity_sliding));

% ============ Calculate AIS Parameters ============
[ais_params] = calculateAISParameters(norm_intensity, PIXEL_TO_MICRON, INTENSITY_FRACTION);

% ============ Plot Results ============
% Create distance array for plotting
distance_array = (1:length(x_pix)) * PIXEL_TO_MICRON;
plotResults(ch1_pic, x_pix, y_pix, intensity_raw, norm_intensity, ...
           distance_array, ais_params, INTENSITY_FRACTION);

% ============ Output Results ============
disp(' ')
disp('AIS Start'), display(ais_params.start)
disp('AIS End'), display(ais_params.end)
disp('AIS Mid'), display(ais_params.mid)
disp('AIS Max'), display(ais_params.max)
disp('AIS Length'), display(ais_params.length)
disp(' ')

clipboard('copy', num2str(ais_params.length, 6));

ais_params.xy_spline = xy_spline;
out = ais_params;

end

% ============ Helper Functions ============

function [x_trace, y_trace] = traceSkeleton(skeleton, start_x, start_y, max_length, max_iterations)
    % Trace skeleton handling branches by finding longest path

    trace_attempts = {};
    iteration = 0;
    branch_found = true;

    while branch_found && iteration < max_iterations
        iteration = iteration + 1;

        % Initialize trace
        x_trace = start_x;
        y_trace = start_y;
        skeleton_copy = skeleton;
        skeleton_copy(start_x, start_y) = 0;

        % Trace until no more pixels or branch found
        branch_found = false;
        while any(skeleton_copy(:))
            % Get current position
            curr_x = x_trace(end);
            curr_y = y_trace(end);

            % Find neighbors
            neighborhood = zeros(size(skeleton_copy));
            x_range = max(1, curr_x-1):min(size(skeleton_copy,1), curr_x+1);
            y_range = max(1, curr_y-1):min(size(skeleton_copy,2), curr_y+1);
            neighborhood(x_range, y_range) = skeleton_copy(x_range, y_range);

            neighbor_idx = find(neighborhood);

            if length(x_trace) > max_length
                break;
            elseif length(neighbor_idx) > 1
                % Branch found - take first neighbor for now
                branch_found = true;
                next_idx = neighbor_idx(1);
                [next_x, next_y] = ind2sub(size(skeleton_copy), next_idx);

                % Mark branch point for later removal
                branch_x = next_x;
                branch_y = next_y;
            elseif length(neighbor_idx) == 1
                [next_x, next_y] = ind2sub(size(skeleton_copy), neighbor_idx);
            else
                break; % No neighbors
            end

            % Add to trace
            x_trace(end+1) = next_x;
            y_trace(end+1) = next_y;
            skeleton_copy(next_x, next_y) = 0;
        end

        % Store trace attempt
        trace_attempts{iteration} = {x_trace, y_trace};

        % If branch found, remove it from skeleton
        if branch_found && exist('branch_x', 'var')
            % Find and remove the branch
            branch_start = find(x_trace == branch_x & y_trace == branch_y, 1);
            if ~isempty(branch_start)
                for i = branch_start:length(x_trace)
                    skeleton(x_trace(i), y_trace(i)) = 0;
                end
            end
            clear branch_x branch_y
        end
    end

    % Select longest trace
    trace_lengths = cellfun(@(x) length(x{1}), trace_attempts);
    [~, longest_idx] = max(trace_lengths);
    x_trace = trace_attempts{longest_idx}{1};
    y_trace = trace_attempts{longest_idx}{2};
end

function [x_pixels, y_pixels, dist_along_axis] = extractUniquePixels(xy_spline, cumulative_dist, pixel_to_micron)
    % Extract unique pixel coordinates along spline
    % CRITICAL: In the original code, the spline is built as [Yais; Xais]
    % So xy_spline(1,:) contains Y coordinates and xy_spline(2,:) contains X coordinates

    % Round spline to whole-number pixel coordinates (matching original)
    xys = round(xy_spline);

    % Extract rounded coordinates with the confusing swap from original
    xs = xys(1,:);  % First row (originally Yais) becomes xs
    ys = xys(2,:);  % Second row (originally Xais) becomes ys

    % Find unique pixels using the exact same logic as original
    double_id = [];
    for i = 1:length(xs)
        for j = (i+1):length(xs)
            if xs(i) == xs(j) && ys(i) == ys(j)
                double_id = [double_id; j];
            end
        end
    end
    double_id = unique(double_id);
    alln = 1:length(xs);
    m = setdiff(alln, double_id);

    % Extract unique coordinates - matching original variable names
    x_pix = xs(m);  % This matches the original x_pix
    y_pix = ys(m);  % This matches the original y_pix

    % Return with the expected names
    x_pixels = x_pix;
    y_pixels = y_pix;

    % Calculate distances for unique pixels (matching original exactly)
    for g = 1:length(x_pixels)
        d = [];
        for h = 1:size(xy_spline, 2)
            d = [d; sqrt((x_pixels(g) - xy_spline(1,h))^2 + (y_pixels(g) - xy_spline(2,h))^2)];
        end
        [~, mindi] = min(d);
        dist_along_axis(g) = cumulative_dist(mindi) * pixel_to_micron;
    end
end

function saveCoordinates(cell_name, x_pixels, y_pixels)
    % Save pixel coordinates to text file

    filename = [cell_name '_xy.txt'];
    fid = fopen(filename, 'wt');
    fprintf(fid, '%f\n', length(x_pixels));
    for i = 1:length(x_pixels)
        fprintf(fid, '%f\t %f\n', x_pixels(i), y_pixels(i));
    end
    fclose(fid);
end

function [intensity_raw, intensity_smooth, intensity_sliding] = extractIntensityProfile(img, x_pixels, y_pixels, pixel_to_micron)
    % Extract fluorescence intensity profile along AIS

    n_pixels = length(x_pixels);
    intensity_raw = zeros(1, n_pixels);
    intensity_smooth = zeros(1, n_pixels);

    % Extract raw and 3x3 smoothed intensities
    for i = 1:n_pixels
        % Raw intensity
        intensity_raw(i) = img(y_pixels(i), x_pixels(i));

        % 3x3 neighborhood average
        y_range = max(1, y_pixels(i)-1):min(size(img,1), y_pixels(i)+1);
        x_range = max(1, x_pixels(i)-1):min(size(img,2), x_pixels(i)+1);
        neighborhood = img(y_range, x_range);
        intensity_smooth(i) = mean(neighborhood(:));
    end

    % Apply sliding window smoothing
    window_pixels = round(1.5 / pixel_to_micron) + 1;
    intensity_sliding = zeros(1, n_pixels);

    for i = 1:n_pixels
        if i < (window_pixels + 1)
            % Start of profile
            window_idx = [1:i, i:min(i+window_pixels, n_pixels)];
        elseif i > (n_pixels - window_pixels - 1)
            % End of profile
            window_idx = [max(1, i-window_pixels):i, i:n_pixels];
        else
            % Middle of profile
            window_idx = (i-window_pixels):(i+window_pixels);
        end
        intensity_sliding(i) = mean(intensity_smooth(window_idx));
    end
end

function ais_params = calculateAISParameters(norm_intensity, pixel_to_micron, intensity_fraction)
    % Calculate AIS start, end, length, and other parameters

    % Find maximum position
    [~, max_idx] = max(norm_intensity);

    % Find start (before max where intensity drops below fraction)
    start_candidates = find((1:length(norm_intensity)) < max_idx & norm_intensity < intensity_fraction);
    if ~isempty(start_candidates)
        start_idx = start_candidates(end);
    else
        start_idx = 1;
    end

    % Find end (after max where intensity drops below fraction)
    end_candidates = find((1:length(norm_intensity)) > max_idx & norm_intensity > intensity_fraction);
    if ~isempty(end_candidates)
        end_idx = end_candidates(end);
    else
        end_idx = length(norm_intensity);
    end

    % Calculate parameters in microns
    ais_params.start = start_idx * pixel_to_micron;
    ais_params.end = end_idx * pixel_to_micron;
    ais_params.length = ais_params.end - ais_params.start;
    ais_params.mid = mean([ais_params.start, ais_params.end]);
    ais_params.max = max_idx * pixel_to_micron;

    % Store indices for plotting
    ais_params.start_idx = start_idx;
    ais_params.end_idx = end_idx;
    ais_params.max_idx = max_idx;
end

function plotResults(img, x_pixels, y_pixels, intensity_raw, norm_intensity, distance_um, ais_params, intensity_fraction)
    % Plot analysis results

    figure(3)

    % Original image with traced path
    subplot(2,2,1)
    imagesc(img)
    colormap(gray)
    axis square
    title('Ch1', 'color', 'g')
    hold on
    plot(x_pixels, y_pixels, '-b')  % Plot with x,y order for display
    hold off

    % Empty subplot for compatibility with original
    subplot(2,2,2)
    axis off

    % Raw intensity profile
    subplot(2,2,3)
    plot(distance_um, intensity_raw, 'g-')
    axis square
    title('Raw')
    xlabel('Distance (μm)')
    ylabel('Intensity')

    % Normalized intensity profile with AIS markers
    subplot(2,2,4)
    plot(distance_um, norm_intensity, 'g-', 'LineWidth', 2)
    hold on
    yline(intensity_fraction, '-.', 'Threshold')

    % Plot AIS boundaries
    plot([ais_params.start ais_params.start], [0 1], 'b-', 'LineWidth', 2)
    plot([ais_params.end ais_params.end], [0 1], 'b-', 'LineWidth', 2)
    plot([ais_params.max ais_params.max], [0 1], 'b-', 'LineWidth', 2)

    text(max(distance_um)-25, 0.9, 'Ch1 prof', 'color', 'b')

    axis square
    title('Smoothed & Normalized')
    xlabel('Distance (μm)')
    ylabel('Normalized Intensity')
    hold off
end