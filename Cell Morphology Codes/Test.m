% Step 1: Load image
img = imread('1_20x.jpg');
gray = rgb2gray(img); % If RGB
bw = imbinarize(gray); % Threshold to binary
% Optional: clean up
bw = imfill(bw, 'holes');
bw = bwareaopen(bw, 50); % Remove small noise
% Step 2: Label connected components
labeled = bwlabel(bw);
% Step 3: Measure properties
props = regionprops(labeled, 'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity');

eccentricity = [props.Eccentricity];
valid_idx = eccentricity > 0.75;
props = props(valid_idx);

% Step 4: Display
imshow(img); 
hold on;
for k = 1:length(props)
    % Get ellipse parameters
    xCenter = props(k).Centroid(1);
    yCenter = props(k).Centroid(2);
    majorAxis = props(k).MajorAxisLength;
    minorAxis = props(k).MinorAxisLength;
    angle = deg2rad(-props(k).Orientation); % Negative to match display rotation
    % Generate ellipse
    t = linspace(0, 2*pi, 100);
    x = (majorAxis/2) * cos(t);
    y = (minorAxis/2) * sin(t);
    R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    coords = R * [x; y];
    xFit = coords(1,:) + xCenter;
    yFit = coords(2,:) + yCenter;
    % Plot ellipse
    plot(xFit, yFit, 'r', 'LineWidth', 1.5);
    % Optional: Print properties
    fprintf('Cell %d - Eccentricity: %.2f\n', k, props(k).Eccentricity);
end
title('Identified Elliptical Cells');