% Step 1: Load and preprocess image
img = imread('1_20x.jpg');
gray = rgb2gray(img);                  % Convert to grayscale
bw = imbinarize(gray);                 % Threshold to binary
bw = imfill(bw, 'holes');              % Fill holes in cells
bw = bwareaopen(bw, 50);               % Remove small noise
% Step 2: Label connected components
labeled = bwlabel(bw);
stats = regionprops(labeled, 'Eccentricity', 'MajorAxisLength', ...
                               'MinorAxisLength', 'Area', 'PixelIdxList');
% Step 3: Create a new binary mask (non-circular only)
eccentricityThreshold = 0.7;   % You can tune this threshold!
filteredMask = false(size(bw));
for i = 1:length(stats)
    ecc = stats(i).Eccentricity;
    if ecc > eccentricityThreshold
        % Keep elliptical (non-circular) cells
        filteredMask(stats(i).PixelIdxList) = true;
    end
end
% Step 4: Visualize result
figure;
subplot(1,3,1); imshow(img); title('Original Image');
subplot(1,3,2); imshow(bw); title('All Cells (Binary)');
subplot(1,3,3); imshow(filteredMask); title('Elliptical Cells Only');
% Step 5: (Optional) Save the output image
imwrite(filteredMask, 'elliptical_cells_only.png');
fprintf('Removed round cells. Retained %d elliptical shapes.\n', bwconncomp(filteredMask).NumObjects);

Detect_Cells_Watershed(gray);

function Detect_Cells_Watershed(gray)
    % Step 1: Normalize image
    norm_image = mat2gray(gray);

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