clear, clc, close all

% ------------------ Step 1: Configuration & Paths ------------------
main_folder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio';
outDir = fullfile(main_folder, 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Parameters
h_max = 12;                 
aspect_ratio_cutoff = 1.2;  
image_cell_diameter = 40;   

% ------------------ Step 2: Filter Files ------------------
all_tifs = dir(fullfile(main_folder, '*.tif'));
pattern = '^1_20_[1-5]\.tif$';
keep = ~cellfun('isempty', regexp({all_tifs.name}, pattern, 'once'));
tif_files = all_tifs(keep);

if isempty(tif_files), error('No matching files found.'); end
cp = cellpose(model="cyto2");

% ------------------ Step 3: Process Images ------------------
for k = 1:length(tif_files)
    fprintf('Processing & Visualizing: %s\n', tif_files(k).name);
    
    img_path = fullfile(tif_files(k).folder, tif_files(k).name);
    img_norm = mean(single(imread(img_path)), 3) ./ 65535;
    img_adj = adapthisteq(img_norm);

    % Segment
    cellLabel = segmentCells2D(cp, img_adj, ImageCellDiameter=image_cell_diameter);
    props = regionprops(cellLabel, img_norm, 'Area', 'MajorAxisLength', 'MinorAxisLength', 'PixelIdxList', 'Centroid');
    
    thicknessMap = zeros(size(img_norm));
    I_bg = median(img_norm(cellLabel == 0));
    
    % Classification and Height Logic
    for i = 1:max(cellLabel(:))
        AR = props(i).MajorAxisLength / props(i).MinorAxisLength;
        idx = props(i).PixelIdxList;
        
        if AR > aspect_ratio_cutoff
            props(i).Class = 'Alive';
            I_cell = max(img_norm(idx));
            cell_h = h_max .* (img_norm(idx) - I_bg) ./ (I_cell - I_bg);
            cell_h(cell_h < 0) = 0;
            thicknessMap(idx) = cell_h;
        else
            props(i).Class = 'Dead';
        end
    end

    % ------------------ Step 4: VISUALS ------------------
    
    % --- Visual 1: 2D Classification Map ---
    % Show original image with color-coded classification
    fig2D = figure('Name', ['Classification - ', tif_files(k).name]);
    imshow(img_adj, []); hold on;
    
    for i = 1:length(props)
        % Get boundary of cell
        mask = (cellLabel == i);
        B = bwboundaries(mask);
        boundary = B{1};
        
        if strcmp(props(i).Class, 'Alive')
            plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1.5); % Green for Alive
        else
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5); % Red for Dead
        end
    end
    title(['Green: Alive | Red: Dead - ', tif_files(k).name]);
    saveas(fig2D, fullfile(outDir, [tif_files(k).name(1:end-4), '_classification.png']));

    % --- Visual 2: 3D Height Map ---
    pbaspect([1 1 0.5]); % take th
    fig3D = figure('Name', ['3D Height - ', tif_files(k).name]);
    [X, Y] = meshgrid(1:size(img_norm,2), 1:size(img_norm,1));
    
    % Use 'surf' but filter out the zeros for cleaner display
    s = surf(X, Y, thicknessMap, 'EdgeColor', 'none');
    colormap(turbo);
    cb = colorbar;
    cb.Label.String = 'Thickness (\mum)';
    
    % Lighting and aesthetic tweaks
    view(45, 35);
    axis tight; grid off;
    camlight headlight; 
    lighting gouraud;
    material shiny;
    
    title(['3D Reconstructed Alive Cells: ', tif_files(k).name]);
    zlabel('Height (\mum)');
    saveas(fig3D, fullfile(outDir, [tif_files(k).name(1:end-4), '_3D_Height.png']));
    
end

fprintf('All visuals generated and saved in: %s\n', outDir);