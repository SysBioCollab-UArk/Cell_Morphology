% this is where I am going to input filenames to get images
filename = '1_20x.jpg';

% function calls below
cell_image = Import_Image(filename);
cell_image = single(cell_image);
cell_image_gray = Convert_to_Grayscale(cell_image);

%Detect_Cell_Morphology(cell_image_gray);
Detect_Cells_Watershed(cell_image_gray);

% functions have to go AFTER function calls...remember that
function cell_image = Import_Image(filename)
    cell_image = imread(filename);
    %imshow(cell_image)
end

% now we are converting to 8-bit grayscale
function cell_image_gray = Convert_to_Grayscale(cell_image)
    cell_image_gray = mean(cell_image, 3);
    %imshow(cell_image_gray)
end

function Detect_Cells_Watershed(cell_image_gray)
    % Step 1: Normalize image
    norm_image = mat2gray(cell_image_gray);

    % Step 2: Invert and smooth
    inverted = imcomplement(imgaussfilt(norm_image, 1.5));  % Try σ = 1.5 for better preservation

    % Step 3: Adaptive thresholding
    BW = imbinarize(inverted, 'adaptive', 'Sensitivity', 0.45);  % Was 0.35
    BW = bwareaopen(BW, 300);     % Was 500: allows smaller cells
    BW = imfill(BW, 'holes');

    % Step 4: Distance transform
    D = -bwdist(~BW);
    D(~BW) = -Inf;

    % Step 5: Marker-controlled watershed
    marker_mask = imextendedmin(D, 1.5);   % Was 2: more sensitive
    D2 = imimposemin(D, marker_mask);
    L = watershed(D2);

    % Step 6: Remove watershed ridges
    BW_watershed = BW;
    BW_watershed(L == 0) = 0;

    % Step 7: Region properties
    stats = regionprops(BW_watershed, 'Centroid', 'Area', ...
        'MajorAxisLength', 'MinorAxisLength');

    % % Step 7: Region properties
    % stats = regionprops(BW_watershed, 'Centroid', 'Eccentricity', 'Area', ...
    %      'MajorAxisLength', 'MinorAxisLength');

    aspect_ratios = [stats.MajorAxisLength] ./ [stats.MinorAxisLength];
    areas = [stats.Area];
    % eccentricity = [stats.Eccentricity];

    % Step 8: Loosen filtering
    valid_idx = areas > 250 & aspect_ratios < 4;  % Looser bounds
    stats = stats(valid_idx);

    % Step 9: Visual overlay
    imshow(labeloverlay(norm_image, BW_watershed));
    hold on;
    for k = 1:length(stats)
        c = stats(k).Centroid;
        text(c(1), c(2), sprintf('%d', k), 'Color', 'y', ...
             'FontSize', 10, 'FontWeight', 'bold');
    end
    hold off;

    % Step 10: Print summary
    fprintf('Detected %d valid cells.\n', length(stats));
    fprintf('Cell # | Area | Aspect Ratio\n');
    for k = 1:length(stats)
        fprintf('%6d | %5.0f | %0.2f\n', k, stats(k).Area, ...
            stats(k).MajorAxisLength / stats(k).MinorAxisLength);
    end
end
