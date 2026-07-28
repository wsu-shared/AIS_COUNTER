function [out] = ais_debug_comparison(cell, threshold, smooth, chosen)
% AIS_DEBUG_COMPARISON - Debug version that saves data for comparison
%
% This version saves intermediate data files that can be compared
% with the original method to find where differences occur

close all

% ============ Constants FROM ORIGINAL ============
PIXEL_TO_MICRON = 0.161;
INTENSITY_FRACTION = 0.33;
MIN_COMPONENT_SIZE = 70;
SMOOTH_KERNEL_SIZE = 20;
SMOOTH_KERNEL_SD = 2;
MORPH_ELEMENT_SIZE = 3;
MAX_AIS_LENGTH = 300;
MAX_BRANCH_ITERATIONS = 10;
SPLINE_SMOOTHING = 0.3;
SLIDING_WINDOW_WIDTH = 1.5;

%% Load images
[filepath, basename, ~] = fileparts(cell);
disp(['Processing: ' cell]);

% Load images exactly as original
filename = [cell '.tif - Processed method 2.5.tif'];
img = imread(filename);
ch1 = img(:,:,1);
ch1_processed = ch1;

% Load original for intensity
ch1_pic = imread([cell '.tif']);

% Display
figure(1)
imagesc(ch1)
colormap(gray)
axis square
hold on
title('Click near AIS start point')

%% Preprocessing (EXACT from original)
img_gray = mat2gray(ch1);

if strcmp(smooth, 'gaussian')
    img_smooth = imfilter(img_gray, fspecial('gaussian', [SMOOTH_KERNEL_SIZE SMOOTH_KERNEL_SIZE], SMOOTH_KERNEL_SD));
else
    img_smooth = img_gray;
end

if threshold == 0
    threshold = graythresh(img_smooth);
end

%% Get click position
disp('Click near AIS start point...')
[x_click, y_click] = ginput(1);
x_click = round(x_click);
y_click = round(y_click);

% Save click position for comparison
save([cell '_debug_click.mat'], 'x_click', 'y_click');
fprintf('Click position: (%d, %d)\n', x_click, y_click);

%% Re-thresholding loop (EXACT from original)
region_selected = false;
loop_count = 0;

while ~region_selected
    % Threshold image
    img_binary = im2bw(img_smooth, threshold);
    
    % Morphological operations
    se_open = strel('square', MORPH_ELEMENT_SIZE);
    se_close = strel('square', MORPH_ELEMENT_SIZE);
    img_opened = imopen(img_binary, se_open);
    img_closed = imclose(img_opened, se_close);
    
    % Find and filter connected components
    cc = bwconncomp(img_closed);
    pixel_counts = cellfun(@numel, cc.PixelIdxList);
    small_components = find(pixel_counts < MIN_COMPONENT_SIZE);
    cc.PixelIdxList(small_components) = [];
    cc.NumObjects = length(cc.PixelIdxList);
    
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
    
    % Check if maximum pixel is in selected region
    if loop_count == 0
        img_double = im2double(ch1);
        max_global = max(img_double(:));
        max_ais = max(max(ais_mask .* img_double));
        
        if max_global == max_ais
            region_selected = true;
            fprintf('Threshold accepted: %.4f\n', threshold);
        else
            disp('Max pixel not in AIS region. Rethresholding...')
            threshold = (max_ais / max_global) * threshold;
            fprintf('New threshold: %.4f\n', threshold);
            loop_count = 1;
        end
    else
        region_selected = true;
    end
end

% Save mask for comparison
save([cell '_debug_mask.mat'], 'ais_mask', 'threshold');

%% Skeletonization
ais_skeleton = bwmorph(ais_mask, 'thin', Inf);

% Save skeleton
save([cell '_debug_skeleton.mat'], 'ais_skeleton');

%% Find starting point
[~, idx_dist] = bwdist(ais_skeleton);
closest_pixel = idx_dist(y_click, x_click);
[start_x, start_y] = ind2sub(size(idx_dist), closest_pixel);

fprintf('Start point: (%d, %d)\n', start_x, start_y);
save([cell '_debug_start.mat'], 'start_x', 'start_y');

%% Trace skeleton
[x_trace, y_trace] = traceSkeleton(ais_skeleton, start_x, start_y, MAX_AIS_LENGTH, MAX_BRANCH_ITERATIONS);

fprintf('Trace length: %d pixels\n', length(x_trace));
save([cell '_debug_trace.mat'], 'x_trace', 'y_trace');

% Display traced path
figure(2)
imagesc(ch1_processed)
colormap(gray)
axis square
title('Traced AIS Path')
hold on
plot(y_trace, x_trace, '.r')

%% Fit spline
xy = double([y_trace; x_trace]);
t = 1:length(xy);
ts = [1:0.1:length(xy), t(end)];
xy_spline = csaps(t, xy, SPLINE_SMOOTHING, ts);

% Calculate distances
distances = zeros(size(xy_spline, 2), 1);
for i = 2:length(xy_spline)
    distances(i) = sqrt((xy_spline(1,i) - xy_spline(1,i-1))^2 + ...
                       (xy_spline(2,i) - xy_spline(2,i-1))^2);
end
cumulative_dist = cumsum(distances);

save([cell '_debug_spline.mat'], 'xy_spline', 'cumulative_dist');

%% Extract unique pixels
[x_pixels, y_pixels, dist_along_axis] = extractUniquePixels(xy_spline, cumulative_dist, PIXEL_TO_MICRON);

fprintf('Unique pixels: %d\n', length(x_pixels));
save([cell '_debug_pixels.mat'], 'x_pixels', 'y_pixels');

% Save coordinates file (like original)
saveCoordinates(cell, x_pixels, y_pixels);

%% Extract intensity profile
[intensity_raw, intensity_smooth, intensity_sliding] = extractIntensityProfile(ch1_pic, x_pixels, y_pixels, PIXEL_TO_MICRON);

% Normalize
norm_intensity = (intensity_sliding - min(intensity_sliding)) / (max(intensity_sliding) - min(intensity_sliding));

save([cell '_debug_intensity.mat'], 'intensity_raw', 'intensity_smooth', ...
     'intensity_sliding', 'norm_intensity');

%% Calculate AIS parameters
ais_params = calculateAISParameters(norm_intensity, PIXEL_TO_MICRON, INTENSITY_FRACTION);

fprintf('\n=== RESULTS ===\n');
fprintf('AIS Start: %.1f µm (idx: %d)\n', ais_params.start, ais_params.start_idx);
fprintf('AIS End: %.1f µm (idx: %d)\n', ais_params.end, ais_params.end_idx);
fprintf('AIS Max: %.1f µm (idx: %d)\n', ais_params.max, ais_params.max_idx);
fprintf('AIS Length: %.1f µm\n', ais_params.length);
fprintf('===============\n');

save([cell '_debug_params.mat'], 'ais_params');

%% Plot results
figure(3)
subplot(2,2,3)
plot(1:length(intensity_raw), intensity_raw, 'g-')
title('Raw Intensity')

subplot(2,2,4)
plot(1:length(norm_intensity), norm_intensity, 'g-', 'LineWidth', 2)
hold on
yline(INTENSITY_FRACTION, '-.', 'Threshold')
plot([ais_params.start_idx ais_params.start_idx], [0 1], 'b-', 'LineWidth', 2)
plot([ais_params.end_idx ais_params.end_idx], [0 1], 'b-', 'LineWidth', 2)
plot([ais_params.max_idx ais_params.max_idx], [0 1], 'r-', 'LineWidth', 2)
title(sprintf('Length = %.1f µm', ais_params.length))

%% Validation click
disp('Click to accept (left) or reject (right/n)...')
[~, ~, button] = ginput(1);

if button == 110 || button == 3  % 'n' or right click
    out = 0;
else
    out = ais_params;
end

% Clean up debug files if accepted
if isstruct(out)
    delete([cell '_debug_*.mat']);
    disp('Debug files cleaned up. Results accepted.');
else
    disp('Results rejected. Debug files kept for analysis.');
end

end

% ============ Helper Functions (EXACT from original) ============

function [x_trace, y_trace] = traceSkeleton(skeleton, start_x, start_y, max_length, max_iterations)
    trace_attempts = {};
    iteration = 0;
    branch_found = true;

    while branch_found && iteration < max_iterations
        iteration = iteration + 1;

        x_trace = start_x;
        y_trace = start_y;
        skeleton_copy = skeleton;
        skeleton_copy(start_x, start_y) = 0;

        branch_found = false;
        while any(skeleton_copy(:))
            curr_x = x_trace(end);
            curr_y = y_trace(end);

            neighborhood = zeros(size(skeleton_copy));
            x_range = max(1, curr_x-1):min(size(skeleton_copy,1), curr_x+1);
            y_range = max(1, curr_y-1):min(size(skeleton_copy,2), curr_y+1);
            neighborhood(x_range, y_range) = skeleton_copy(x_range, y_range);

            neighbor_idx = find(neighborhood);

            if length(x_trace) > max_length
                break;
            elseif length(neighbor_idx) > 1
                branch_found = true;
                next_idx = neighbor_idx(1);
                [next_x, next_y] = ind2sub(size(skeleton_copy), next_idx);
                branch_x = next_x;
                branch_y = next_y;
            elseif length(neighbor_idx) == 1
                [next_x, next_y] = ind2sub(size(skeleton_copy), neighbor_idx);
            else
                break;
            end

            x_trace(end+1) = next_x;
            y_trace(end+1) = next_y;
            skeleton_copy(next_x, next_y) = 0;
        end

        trace_attempts{iteration} = {x_trace, y_trace};

        if branch_found && exist('branch_x', 'var')
            branch_start = find(x_trace == branch_x & y_trace == branch_y, 1);
            if ~isempty(branch_start)
                for i = branch_start:length(x_trace)
                    skeleton(x_trace(i), y_trace(i)) = 0;
                end
            end
            clear branch_x branch_y
        end
    end

    trace_lengths = cellfun(@(x) length(x{1}), trace_attempts);
    [~, longest_idx] = max(trace_lengths);
    x_trace = trace_attempts{longest_idx}{1};
    y_trace = trace_attempts{longest_idx}{2};
end

function [x_pixels, y_pixels, dist_along_axis] = extractUniquePixels(xy_spline, cumulative_dist, pixel_to_micron)
    xys = round(xy_spline);
    xs = xys(1,:);
    ys = xys(2,:);

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

    x_pix = xs(m);
    y_pix = ys(m);

    x_pixels = x_pix;
    y_pixels = y_pix;

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
    filename = [cell_name '_xy.txt'];
    fid = fopen(filename, 'wt');
    fprintf(fid, '%f\n', length(x_pixels));
    for i = 1:length(x_pixels)
        fprintf(fid, '%f\t %f\n', x_pixels(i), y_pixels(i));
    end
    fclose(fid);
end

function [intensity_raw, intensity_smooth, intensity_sliding] = extractIntensityProfile(img, x_pixels, y_pixels, pixel_to_micron)
    n_pixels = length(x_pixels);
    intensity_raw = zeros(1, n_pixels);
    intensity_smooth = zeros(1, n_pixels);

    for i = 1:n_pixels
        intensity_raw(i) = img(y_pixels(i), x_pixels(i));

        y_range = max(1, y_pixels(i)-1):min(size(img,1), y_pixels(i)+1);
        x_range = max(1, x_pixels(i)-1):min(size(img,2), x_pixels(i)+1);
        neighborhood = img(y_range, x_range);
        intensity_smooth(i) = mean(neighborhood(:));
    end

    window_pixels = round(1.5 / pixel_to_micron) + 1;
    intensity_sliding = zeros(1, n_pixels);

    for i = 1:n_pixels
        if i < (window_pixels + 1)
            window_idx = [1:i, i:min(i+window_pixels, n_pixels)];
        elseif i > (n_pixels - window_pixels - 1)
            window_idx = [max(1, i-window_pixels):i, i:n_pixels];
        else
            window_idx = (i-window_pixels):(i+window_pixels);
        end
        intensity_sliding(i) = mean(intensity_smooth(window_idx));
    end
end

function ais_params = calculateAISParameters(norm_intensity, pixel_to_micron, intensity_fraction)
    [~, max_idx] = max(norm_intensity);

    start_candidates = find((1:length(norm_intensity)) < max_idx & norm_intensity < intensity_fraction);
    if ~isempty(start_candidates)
        start_idx = start_candidates(end);
    else
        start_idx = 1;
    end

    end_candidates = find((1:length(norm_intensity)) > max_idx & norm_intensity > intensity_fraction);
    if ~isempty(end_candidates)
        end_idx = end_candidates(end);
    else
        end_idx = length(norm_intensity);
    end

    ais_params.start = start_idx * pixel_to_micron;
    ais_params.end = end_idx * pixel_to_micron;
    ais_params.length = ais_params.end - ais_params.start;
    ais_params.mid = mean([ais_params.start, ais_params.end]);
    ais_params.max = max_idx * pixel_to_micron;

    ais_params.start_idx = start_idx;
    ais_params.end_idx = end_idx;
    ais_params.max_idx = max_idx;
end