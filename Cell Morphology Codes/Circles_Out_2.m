% Load and preprocess the image
img = imread('1_20x.jpg');

if size(img, 3) == 3
    gray = rgb2gray(img);
else
    gray = img;
end

Detect_Elliptical_Cells(gray, img);

function Detect_Elliptical_Cells(gray, originalImage)
    % Step 1: Normalize grayscale image
    norm_image = mat2gray(gray);

    % Step 2: Invert and smooth (less aggressive blur)
    inverted = imcomplement(imgaussfilt(norm_image, 1.0));  % σ = 1.0

    % Step 3: Adaptive thresholding
    BW = imbinarize(inverted, 'adaptive', 'Sensitivity', 0.5);  % was 0.45
    BW = bwareaopen(BW, 300);     % Allow small objects
    BW = imfill(BW, 'holes');
    BW = imdilate(BW, strel('disk', 1));  % Fill thin parts of elongated cells

    % Step 4: Distance transform
    D = -bwdist(~BW);
    D(~BW) = -Inf;

    % Step 5: Marker-controlled watershed
    marker_mask = imextendedmin(D, 1.5);   % Sensitive minima
    D2 = imimposemin(D, marker_mask);
    L = watershed(D2);

    % Step 6: Remove watershed ridges
    BW_watershed = BW;
    BW_watershed(L == 0) = 0;

    % Step 7: Region properties
    stats = regionprops(BW_watershed, 'Centroid', 'Area', ...
        'MajorAxisLength', 'MinorAxisLength', 'Eccentricity');

    % Extract features
    areas = [stats.Area];
    aspect_ratios = [stats.MajorAxisLength] ./ [stats.MinorAxisLength];
    eccentricity = [stats.Eccentricity];

    % Step 8: Filter elliptical cells
    valid_idx = (areas > 250) & (aspect_ratios < 4) & (eccentricity > 0.75);
    stats = stats(valid_idx);

    % Step 9: Visual overlay
    figure;
    imshow(labeloverlay(originalImage, BW_watershed));
    hold on;
    for k = 1:length(stats)
        c = stats(k).Centroid;
        text(c(1), c(2), sprintf('%d', k), 'Color', 'y', ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
    hold off;
    title('Detected Elliptical Cells');

    % Step 10: Print summary
    fprintf('Detected %d valid elliptical cells.\n', length(stats));
    fprintf('Cell # | Area | Aspect Ratio | Eccentricity\n');
    for k = 1:length(stats)
        fprintf('%6d | %5.0f | %0.2f | %0.2f\n', k, stats(k).Area, ...
            stats(k).MajorAxisLength / stats(k).MinorAxisLength, ...
            stats(k).Eccentricity);
    end

    % Step 11: Save result
    imwrite(BW_watershed, 'elliptical_cells_only.png');
end