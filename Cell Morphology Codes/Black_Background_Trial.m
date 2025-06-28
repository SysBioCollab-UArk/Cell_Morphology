% Load grayscale image
img = imread('1_20x_W.tif');
if size(img,3) == 3
    gray = rgb2gray(img);
else
    gray = img;
end

% Step 1: Enhance contrast
gray_eq = adapthisteq(gray);

% Step 2: Gradient magnitude to find outlines
[Gmag, ~] = imgradient(gray_eq);
Gmag = mat2gray(Gmag);  % Normalize to [0, 1]

% Step 3: Threshold strong gradient regions (halos)
halo_mask = imbinarize(Gmag, 'adaptive', ...
    'ForegroundPolarity', 'bright', ...
    'Sensitivity', 0.45);
halo_mask = imclose(halo_mask, strel('disk', 2));
halo_mask = imfill(halo_mask, 'holes');
halo_mask = bwareaopen(halo_mask, 100);
figure; imshow(halo_mask); title('Gradient-Based Halo Candidates');

% Step 4: Label and filter for round objects
labeled = bwlabel(halo_mask);
stats = regionprops(labeled, 'Area', 'Eccentricity', ...
                    'Solidity', 'MajorAxisLength', 'MinorAxisLength');

figure;
imshow(label2rgb(labeled));
title('All Detected Regions (Before Filtering)');

% Keep: round, solid, medium-area blobs
% DEBUG: How many candidates before filtering?
fprintf('Total objects detected: %d\n', length(stats));

% Loosen the filtering temporarily to avoid excluding everything
valid = find([stats.Area] > 100);  % Keep anything moderately sized

% DEBUG: How many remain after filtering?
fprintf('Objects kept after filtering: %d\n', length(valid));

% Use mask only if something was kept
if ~isempty(valid)
    cell_mask = ismember(labeled, valid);
else
    warning('No valid regions detected. Falling back to unfiltered mask.');
    cell_mask = halo_mask;  % Use full mask before filtering
end

cell_mask = ismember(labeled, valid);
figure; imshow(cell_mask); title('Filtered Cell Mask (Halo Shape + Solidity)');

% Step 5: Mask the original image
masked_img = img;
if nnz(cell_mask) > 0
    if size(img,3) == 3
        for c = 1:3
            ch = masked_img(:,:,c);
            ch(~cell_mask) = 0;
            masked_img(:,:,c) = ch;
        end
    else
        masked_img(~cell_mask) = 0;
    end
end

figure;
imshow(masked_img);
title('Final Masked Cells (Morphology Filtered)');