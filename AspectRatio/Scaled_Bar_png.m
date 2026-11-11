
%% Batch Add Scale Bars to Microscopy Images

clear;
clc;
close all;

% ---------------------------------------------------------
% Step 1: Define Main Folder
% ---------------------------------------------------------

mainFolder = ...
    '/Users/alexandravega/PycharmProjects/Cell_Morphology/AspectRatio/BioRep5_07_13_2026';

% ---------------------------------------------------------
% Step 2: Create Output Folder
% ---------------------------------------------------------
% Original PNG files will NOT be modified.

outputFolder = fullfile(mainFolder, 'Images_With_ScaleBars');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% ---------------------------------------------------------
% Step 3: Find All PNG Images
% ---------------------------------------------------------

pngFiles = dir(fullfile(mainFolder, '*.png'));

if isempty(pngFiles)
    error('No PNG images were found in: %s', mainFolder);
end

fprintf('\n========================================\n');
fprintf('SCALE BAR BATCH PROCESSING\n');
fprintf('========================================\n');
fprintf('Images found: %d\n', length(pngFiles));
fprintf('Output folder:\n%s\n\n', outputFolder);

% ---------------------------------------------------------
% Step 4: Process Each PNG Image
% ---------------------------------------------------------

for i = 1:length(pngFiles)

    fileName = pngFiles(i).name;
    imagePath = fullfile(mainFolder, fileName);

    fprintf('----------------------------------------\n');
    fprintf('Processing: %s\n', fileName);

    % -----------------------------------------------------
    % Read image
    % -----------------------------------------------------

    img = imread(imagePath);

    % -----------------------------------------------------
    % Read PNG metadata
    % -----------------------------------------------------

    info = imfinfo(imagePath);

    % -----------------------------------------------------
    % Check Calibration Metadata
    % -----------------------------------------------------

    if isempty(info.XResolution) || isempty(info.YResolution)

        warning('No spatial resolution information found for %s. Image skipped.', ...
            fileName);

        continue;

    end

    % -----------------------------------------------------
    % Determine Micrometers per Pixel
    % -----------------------------------------------------

    if strcmpi(info.ResolutionUnit, 'Centimeter')

        % 1 centimeter = 10,000 micrometers

        micronsPerPixelX = 10000 / info.XResolution;
        micronsPerPixelY = 10000 / info.YResolution;

    elseif strcmpi(info.ResolutionUnit, 'Inch')

        % 1 inch = 25,400 micrometers

        micronsPerPixelX = 25400 / info.XResolution;
        micronsPerPixelY = 25400 / info.YResolution;

    elseif strcmpi(info.ResolutionUnit, 'Meter')

        % 1 meter = 1,000,000 micrometers

        micronsPerPixelX = 1000000 / info.XResolution;
        micronsPerPixelY = 1000000 / info.YResolution;

    else

        warning(['Unsupported resolution unit "%s" for %s. ' ...
                 'Image skipped.'], ...
                 info.ResolutionUnit, fileName);

        continue;

    end

    % -----------------------------------------------------
    % Check X and Y Calibration
    % -----------------------------------------------------

    calibrationDifference = ...
        abs(micronsPerPixelX - micronsPerPixelY);

    if calibrationDifference > 0.001

        warning(['X and Y calibration differ for %s. ' ...
                 'Please verify this image.'], ...
                 fileName);

    end

    % Use X calibration for horizontal scale bar

    micronsPerPixel = micronsPerPixelX;

    % -----------------------------------------------------
    % Image Dimensions
    % -----------------------------------------------------

    [imageHeight, imageWidth, ~] = size(img);

    % -----------------------------------------------------
    % Calculate Physical Field of View
    % -----------------------------------------------------

    fieldWidthMicrons = ...
        imageWidth * micronsPerPixelX;

    fieldHeightMicrons = ...
        imageHeight * micronsPerPixelY;

    % -----------------------------------------------------
    % Automatically Select Scale-Bar Length
    % -----------------------------------------------------
    %
    % Target approximately 15% of image width.
    %
    % The program then selects a conventional rounded
    % scale-bar value.

    targetScaleBarMicrons = ...
        0.15 * fieldWidthMicrons;

    niceValues = ...
        [1 2 5 10 20 25 50 100 200 250 500 1000 2000];

    possibleValues = ...
        niceValues(niceValues <= targetScaleBarMicrons);

    if isempty(possibleValues)

        scaleBarMicrons = niceValues(1);

    else

        scaleBarMicrons = possibleValues(end);

    end

    % -----------------------------------------------------
    % Convert Scale-Bar Length from Micrometers to Pixels
    % -----------------------------------------------------

    scaleBarPixels = ...
        scaleBarMicrons / micronsPerPixel;

    % -----------------------------------------------------
    % Scale Position According to Image Dimensions
    % -----------------------------------------------------
    %
    % Using percentages rather than fixed pixel values
    % allows the script to work with images of different
    % sizes.

    marginRight = round(imageWidth * 0.05);
    marginBottom = round(imageHeight * 0.08);

    x2 = imageWidth - marginRight;
    x1 = x2 - scaleBarPixels;

    y = imageHeight - marginBottom;

    % -----------------------------------------------------
    % Tick Height and Label Position
    % -----------------------------------------------------

    tickHeight = round(imageHeight * 0.012);

    labelOffset = round(imageHeight * 0.035);

    % -----------------------------------------------------
    % Create Figure
    % -----------------------------------------------------

    fig = figure( ...
        'Visible', 'off', ...
        'Color', 'white');

    % -----------------------------------------------------
    % Display Image
    % -----------------------------------------------------

    if ndims(img) == 2

        % Grayscale microscopy image
        imshow(img, []);

    else

        % RGB image
        imshow(img);

    end

    hold on;

    % -----------------------------------------------------
    % Draw Horizontal Scale Bar
    % -----------------------------------------------------

    line([x1 x2], [y y], ...
        'Color', 'k', ...
        'LineWidth', 4);

    % -----------------------------------------------------
    % Draw Left Vertical Tick
    % -----------------------------------------------------

    line([x1 x1], ...
        [y-tickHeight y+tickHeight], ...
        'Color', 'k', ...
        'LineWidth', 3);

    % -----------------------------------------------------
    % Draw Right Vertical Tick
    % -----------------------------------------------------

    line([x2 x2], ...
        [y-tickHeight y+tickHeight], ...
        'Color', 'k', ...
        'LineWidth', 3);

    % -----------------------------------------------------
    % Add Scale-Bar Label
    % -----------------------------------------------------

    text((x1+x2)/2, ...
        y-labelOffset, ...
        sprintf('%g \\mum', scaleBarMicrons), ...
        'Color', 'k', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');

    axis image;
    axis off;

    % -----------------------------------------------------
    % Save Annotated Image
    % -----------------------------------------------------

    [~, baseName, ~] = fileparts(fileName);

    outputFile = fullfile( ...
        outputFolder, ...
        [baseName '_ScaleBar.png']);

    exportgraphics(gca, ...
        outputFile, ...
        'Resolution', 300);

    % -----------------------------------------------------
    % Close Figure
    % -----------------------------------------------------

    close(fig);

    % -----------------------------------------------------
    % Display Calibration Information
    % -----------------------------------------------------

    fprintf('Image size: %d x %d pixels\n', ...
        imageWidth, imageHeight);

    fprintf('Calibration X: %.6f um/pixel\n', ...
        micronsPerPixelX);

    fprintf('Calibration Y: %.6f um/pixel\n', ...
        micronsPerPixelY);

    fprintf('Field of view: %.2f x %.2f um\n', ...
        fieldWidthMicrons, ...
        fieldHeightMicrons);

    fprintf('Selected scale bar: %g um\n', ...
        scaleBarMicrons);

    fprintf('Scale bar length: %.2f pixels\n', ...
        scaleBarPixels);

    fprintf('Saved as: %s\n\n', ...
        outputFile);

end

% ---------------------------------------------------------
% Step 5: Finished
% ---------------------------------------------------------

fprintf('========================================\n');
fprintf('BATCH PROCESSING COMPLETE\n');
fprintf('========================================\n');

fprintf('Annotated images were saved in:\n%s\n', ...
    outputFolder);
