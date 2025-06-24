% Step 1: Load image and preprocess
img = imread('1_20x.jpg');
gray = rgb2gray(img);              % Convert to grayscale if RGB
bw = imbinarize(gray);             % Convert to binary using thresholding

% Optional cleanup
bw = imfill(bw, 'holes');          % Fill internal holes in cells
bw = bwareaopen(bw, 50);           % Remove small noise (< 50 px)

% Step 2: Label connected components
labeled = bwlabel(bw);

% Step 3: Measure cell properties
props = regionprops(labeled, 'Centroid', 'Orientation', ...
    'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Area');

% Filter by eccentricity
eccentricity = [props.Eccentricity];
valid_idx = eccentricity > 0.9;   % Keep mostly elongated cells
props = props(valid_idx);          % Update to only valid cells

% Step 4: Display original image with ellipses
figure;                            % Open a new figure window
imshow(gray);                       % Show the original image
hold on;                           % Allow drawing on top
axis image;                        % Keep correct aspect ratio
axis off;                          % Optional: hide axis ticks

for k = 1:length(props)
    % Get ellipse parameters
    xCenter = props(k).Centroid(1);
    yCenter = props(k).Centroid(2);
    majorAxis = props(k).MajorAxisLength;
    minorAxis = props(k).MinorAxisLength;
    angle = deg2rad(-props(k).Orientation);  % Convert orientation to radians

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

    % Print properties
    fprintf('Cell %d - Eccentricity: %.2f\n', k, props(k).Eccentricity);
end

title('Identified Elliptical Cells');

aspect_ratios = [props.MajorAxisLength] ./ [props.MinorAxisLength];

for k = 1:length(props)
        fprintf('%6d | %5.0f | %0.2f\n', k, props(k).Area, ...
            props(k).MajorAxisLength / props(k).MinorAxisLength);
end