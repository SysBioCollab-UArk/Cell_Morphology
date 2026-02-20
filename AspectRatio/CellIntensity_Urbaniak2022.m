clear, clc, close all

% ------------------ Step 1: Set main folder ------------------
main_folder = '/Users/alexandravega/PycharmProjects/Cell_Morphology';
% ------------------ Step 2: Get all TIF images ------------------
tif_files = dir([fullfile(main_folder), '/AspectRatio/*.tif']);

cp = cellpose(model="cyto2");

%% Load and preprocess Image  ( BEFORE)
%img = imread('cell_phase_contrast.tif');
%img =imread([tif_files(1).folder, '/', tif_files(1).name]);

%if size(img,3)==3
    img = rgb2gray(img);
%end

%img = im2double(img);
%img = imadjust(img);

%% Load and preprocess image
img = imread(fullfile(tif_files(1).folder, tif_files(1).name));

% Convert to grayscale if multi-channel
if ndims(img) == 3
    img = mean(single(img), 3);   % safer than rgb2gray for microscopy TIFFs
else
    img = single(img);
end

% Normalize intensity (assumes 16-bit microscopy image)
img = img ./ 65535;         % robust to 8-bit or 16-bit data

% Contrast enhancement (mild, segmentation-friendly)
img_adj = adapthisteq(img);

%% Segment cells 
%{
% Adaptive thresholding works best for phase contrast
bw = imbinarize(img,'adaptive','Sensitivity',0.45);

% Clean segmentation
bw = bwareaopen(bw,200);
bw = imfill(bw,'holes');

figure; imshow(bw); title('Cell mask');

% Label connected components in the binary image
[labeledCells, numCells] = bwlabel(bw);
%}

cellLabel = segmentCells2D(cp, img_adj, ImageCellDiameter=40);

%% Get intensity of each cell and the background intensity
% Get the max intensity of the cells and the median intensity of the background
I_bg = median(img(cellLabel == 0));

cellIntensity = zeros(1,max(cellLabel(:)));

for i = 1:max(cellLabel(:))
    % Find where the cell is
    mask = (cellLabel == i);

    % Calculate average intensity within a cell
    cellIntensity(i) = mean(img(mask));
end

I_cell = max(cellIntensity);

%% Thickness reconstruction
h_max = 12;   % micrometers (example; adjust for your cell type)

% Linear intensity → thickness conversion
thickness = h_max .* ...
    (img - I_bg) / (I_cell - I_bg);

% Ensure no negative values
thickness(thickness < 0) = 0;

%% Volume estimation
cellVolume = zeros(1,max(cellLabel(:)));

for i = 1:max(cellLabel(:))
    % Find where the cell is
    mask = (cellLabel == i);

    % Calculate average intensity within a cell
    cellVolume(i) = sum(thickness(mask));
end

%% Visualization

[X,Y] = meshgrid(1:size(img,2),1:size(img,1));

figure;
surf(X,Y,thickness,'EdgeColor','none');
colormap turbo
view(45,30)
axis tight
xlabel('X (px)'); ylabel('Y (px)'); zlabel('Thickness (µm)');
title('Reconstructed Cell Thickness Map');

%% Extract quantitative thickness metrics

mean_thickness = mean(thickness(bw));
max_thickness  = max(thickness(bw));
cell_volume = sum(thickness(bw));  % in pixel·µm units

fprintf('Mean thickness: %.2f µm\n', mean_thickness);
fprintf('Max thickness: %.2f µm\n', max_thickness);
 
% Save the thickness map as a TIF image
imwrite(thickness, fullfile(main_folder, 'thickness_map.tif'));

