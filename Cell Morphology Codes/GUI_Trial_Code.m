% this is where I am going to input filenames to get images
filename = '1_20x.jpg';

% function calls for function below
cell_image = Import_Image(filename);
cell_image = single(cell_image);
cell_image_gray = Convert_to_Grayscale(cell_image);

%Detect_Cell_Morphology(cell_image_gray);
%Detect_Cells_Watershed(cell_image_gray);
CellDetectionGUI();


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

function CellDetectionGUI()
    % Load and convert image
    filename = '1_20x.jpg';
    cell_image = single(imread(filename));
    cell_image_gray = Convert_to_Grayscale(cell_image);

    % Create figure
    f = figure('Name','Cell Detector GUI','Position',[100 100 1000 600]);

    % Axes for image
    ax = axes('Parent',f,'Position',[0.05 0.2 0.6 0.75]);

    % Sliders for parameters
    uicontrol('Style','text','Position',[700 510 120 20],'String','Sensitivity');
    s1 = uicontrol('Style','slider','Min',0.2,'Max',0.8,'Value',0.45,...
        'Position',[700 490 150 20],'Callback',@updateDisplay);

    uicontrol('Style','text','Position',[700 450 120 20],'String','Min Area');
    s2 = uicontrol('Style','slider','Min',100,'Max',1000,'Value',300,...
        'Position',[700 430 150 20],'Callback',@updateDisplay);

    uicontrol('Style','text','Position',[700 390 120 20],'String','Max Aspect');
    s3 = uicontrol('Style','slider','Min',1.5,'Max',6,'Value',4,...
        'Position',[700 370 150 20],'Callback',@updateDisplay);

    % Run initial detection
    updateDisplay();

    function updateDisplay(~,~)
        sensitivity = s1.Value;
        min_area = round(s2.Value);
        max_aspect = s3.Value;

        stats = Detect_Cells_GUI(cell_image_gray, sensitivity, min_area, max_aspect, ax);

        fprintf('Updated: %d cells\n', numel(stats));
    end
end

function stats = Detect_Cells_GUI(cell_image_gray, sensitivity, area_thresh, aspect_thresh, ax)
    norm_image = mat2gray(cell_image_gray);
    norm_image = adapthisteq(norm_image);
    inverted = imcomplement(imgaussfilt(norm_image, 1.5));
    BW = imbinarize(inverted, 'adaptive', 'Sensitivity', sensitivity);
    BW = bwareaopen(BW, 300);
    BW = imfill(BW, 'holes');
    D = -bwdist(~BW);
    D(~BW) = -Inf;
    marker_mask = imregionalmax(imgaussfilt(D, 2));
    D2 = imimposemin(D, marker_mask);
    L = watershed(D2);
    BW_watershed = BW;
    BW_watershed(L == 0) = 0;
    stats = regionprops(BW_watershed, 'Centroid', 'Area', ...
        'MajorAxisLength', 'MinorAxisLength');
    aspect_ratios = [stats.MajorAxisLength] ./ [stats.MinorAxisLength];
    areas = [stats.Area];
    valid_idx = areas > area_thresh & aspect_ratios < aspect_thresh;
    stats = stats(valid_idx);

    % Display on axes
    axes(ax); cla(ax);
    imshow(labeloverlay(norm_image, BW_watershed), 'Parent', ax); hold on;
    for k = 1:length(stats)
        c = stats(k).Centroid;
        text(c(1), c(2), sprintf('%d', k), 'Color', 'y', ...
             'FontSize', 10, 'FontWeight', 'bold');
    end
    hold off;
end

