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



% Step 1: Label filtered elliptical cells
labeledEllipses = bwlabel(filteredMask);
% Step 2: Get major axis info
props = regionprops(labeledEllipses, 'Centroid', 'Orientation', 'MajorAxisLength');
% Step 3: Display original image and overlay lines
imshow(img); hold on;
title('Elliptical Cells with Major Axis Lines');
for i = 1:length(props)
    % Extract properties
    center = props(i).Centroid;
    angle_deg = props(i).Orientation;
    angle_rad = deg2rad(-angle_deg); % Negative because image Y axis goes down
    majorLen = props(i).MajorAxisLength;
    % Compute endpoints of the major axis line
    dx = (majorLen / 2) * cos(angle_rad);
    dy = (majorLen / 2) * sin(angle_rad);
    x1 = center(1) - dx;
    y1 = center(2) - dy;
    x2 = center(1) + dx;
    y2 = center(2) + dy;
    % Plot the line
    plot([x1, x2], [y1, y2], 'g-', 'LineWidth', 2);
    % (Optional) Mark the center
    plot(center(1), center(2), 'b+', 'MarkerSize', 5, 'LineWidth', 1.5);
end