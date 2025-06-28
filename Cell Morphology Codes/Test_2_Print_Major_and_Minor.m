% Step 1: Load image and preprocess
filename = ('1_20x.jpg');
img = imread(filename);
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

fprintf(' Index |   Area  | MajorAx | MinorAx | Aspect\n');
fprintf('-------|---------|---------|---------|--------\n');

for k = 1:length(props)
    fprintf('%6d | %7.0f | %7.2f | %7.2f | %7.2f\n', ...
        k, ...
        props(k).Area, ...
        props(k).MajorAxisLength, ...
        props(k).MinorAxisLength, ...
        props(k).MajorAxisLength / props(k).MinorAxisLength);
end

% Step 5: Export properties to Excel
T = table((1:length(props))', ...
          [props.Area]', ...
          [props.MajorAxisLength]', ...
          [props.MinorAxisLength]', ...
          [props.MajorAxisLength]' ./ [props.MinorAxisLength]', ...
          'VariableNames', {'Index', 'Area', 'MajorAxis', 'MinorAxis', 'AspectRatio'});

% % Write to Excel file
% writetable(T, 'cell_properties.xlsx');
% fprintf('\nCell properties successfully written to "cell_properties.xlsx".\n');
% 
% winopen('cell_properties.xlsx');

% Derive Excel filename from image filename
[~, baseName, ~] = fileparts(filename);  % Extract base name, no extension
excelName = [baseName '_cell_properties.xlsx'];  % e.g., "1_20x_cell_properties.xlsx"

% Write to Excel file
writetable(T, excelName);
fprintf('\nCell properties successfully written to "%s".\n', excelName);

% Open Excel file (Windows only)
winopen(excelName);
