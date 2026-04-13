clear, clc, close all
% ------------------ Step 1: Configuration & Paths ------------------
main_folder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026';
outDir = fullfile(main_folder, 'output');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Parameters
h_max = 12;                 
aspect_ratio_cutoff = 1.1;  
circ_cutoff = 0.98;         
image_cell_diameter = 40;   
pixelSize = 0.5357154;      
pixelArea = pixelSize^2;    

% ------------------ Step 2: Filter Files ------------------
all_tifs = dir(fullfile(main_folder, '*.tif'));
pattern = '^1_[1-2]\.tif$';
keep = ~cellfun('isempty', regexp({all_tifs.name}, pattern, 'once'));
tif_files = all_tifs(keep);
if isempty(tif_files), error('No matching files found.'); end
cp = cellpose(model="cyto2");

% ------------------ Step 3: Process Images ------------------
for k = 1:length(tif_files)
    fprintf('Processing: %s\n', tif_files(k).name);
    
    img_path = fullfile(tif_files(k).folder, tif_files(k).name);
    img_raw = imread(img_path);
    img_norm = mean(single(img_raw), 3) ./ 65535;
    img_adj = adapthisteq(img_norm);
    
    % --- LOOP FIX: The analysis now happens INSIDE the k-loop ---
    
    % Segment
    cellLabel = segmentCells2D(cp, img_adj, ImageCellDiameter=image_cell_diameter);
    props = regionprops(cellLabel, img_norm, 'Area', 'MajorAxisLength', ...
        'MinorAxisLength', 'PixelIdxList', 'Centroid', 'Orientation', 'Circularity');
    
    % Identify Border Cells
    borderMask = false(size(cellLabel));
    borderMask([1, end], :) = true;
    borderMask(:, [1, end]) = true;
    borderIDs = unique(cellLabel(borderMask));
    borderIDs(borderIDs == 0) = []; 
    
    thicknessMap = zeros(size(img_norm));
    I_bg = median(img_norm(cellLabel == 0));
    csvData = table();
    
    % Classification, Volume, and Table Logic
    for i = 1:length(props)
        aspectRatio = props(i).MajorAxisLength / props(i).MinorAxisLength;
        idx = props(i).PixelIdxList;
        isBorder = ismember(i, borderIDs);
        passMorphology = (aspectRatio > aspect_ratio_cutoff && props(i).Circularity < circ_cutoff);
        
        currentVolume = 0; 
        
        if isBorder
            status = 'Border Exclusion';
            textColor = 'red';
        elseif ~passMorphology
            status = 'Dead/Debris';
            textColor = 'red';
        else
            status = 'Viable Cells';
            textColor = 'green';
            
            I_cell = max(img_norm(idx));
            cell_h = h_max .* (img_norm(idx) - I_bg) ./ (I_cell - I_bg);
            cell_h(cell_h < 0) = 0;
            thicknessMap(idx) = cell_h;
            currentVolume = sum(cell_h) * pixelArea;
        end
        
        props(i).DisplayColor = textColor;
        
        newRow = table(i, props(i).Area, props(i).MajorAxisLength * pixelSize, ...
            props(i).MinorAxisLength * pixelSize, aspectRatio, props(i).Orientation, ...
            props(i).Circularity, currentVolume, props(i).Centroid(1), ...
            props(i).Centroid(2), {status}, 'VariableNames', ...
            {'Cell_Number', 'Area_px', 'MajorAxis_um', 'MinorAxis_um', ...
            'AspectRatio', 'Orientation', 'Circularity', 'Volume_um3', ...
            'Centroid_X', 'Centroid_Y', 'Status'});
        csvData = [csvData; newRow];
    end
    
    % Save CSV
    csv_name = [tif_files(k).name(1:end-4), '_morphology.csv'];
    writetable(csvData, fullfile(outDir, csv_name));

    % ------------------ Step 4: VISUALS ------------------
    
    % Visual 1: 2D Classification Map
    fig2D = figure('Name', ['Classification - ', tif_files(k).name], 'Visible', 'off');
    imshow(img_adj, []); hold on;
    for i = 1:length(props)
        mask = (cellLabel == i);
        B = bwboundaries(mask);
        if ~isempty(B)
            boundary = B{1};
            plot(boundary(:,2), boundary(:,1), props(i).DisplayColor, 'LineWidth', 1.5);
        end
    end
    title(['Green: Viable | Red: Excluded - ', tif_files(k).name]);
    saveas(fig2D, fullfile(outDir, [tif_files(k).name(1:end-4), '_classification.png']));
    close(fig2D); 
    
    % Visual 2: INTERACTIVE 3D Height Map (BLACK BACKGROUND)
    fig3D = figure('Name', ['3D Interactive - ', tif_files(k).name], 'Color', 'k'); % Black figure
    ax = axes('Parent', fig3D, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); % Black axes, white labels
    [X, Y] = meshgrid(1:size(img_norm,2), 1:size(img_norm,1));
    
    surf(ax, X, Y, thicknessMap, 'EdgeColor', 'none');
    colormap(turbo); 
    cb = colorbar;
    cb.Color = 'w'; % White text for colorbar
    cb.Label.String = 'Thickness (\mum)';
    cb.Label.Color = 'w';
    
    view(45, 35); 
    axis tight; 
    grid on; 
    set(ax, 'GridColor', [0.3 0.3 0.3]); % Subtle grey grid lines
    
    camlight headlight; 
    lighting gouraud;
    
    title(['3D Reconstruction: ', tif_files(k).name], 'Color', 'w');
    zlabel('Height (\mum)', 'Color', 'w');
    pbaspect([1 1 0.5]);
    
    rotate3d on; 
end 
fprintf('Done. All images processed.\n');