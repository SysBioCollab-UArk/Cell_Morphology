clear, clc, close all

% ------------------ Step 1: Configuration & Paths ------------------
main_folder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026';
%main_folder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/2026_04_09_10_38_48--BioRep2_Double_Cells_04_09_2026_Final';

outDir = fullfile(main_folder, 'Urbaniak_Volumes');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Parameters
h_max = 12;                       
image_cell_diameter = 40;   
pixelSize = 0.5357154;      
pixelArea = pixelSize^2;    

% ------------------ Step 2: Filter Files ------------------
all_tifs = dir(fullfile(main_folder, '*.tif'));

% Update the pattern to match ANY number_number.tif format (e.g., 1_1, 5_4, 12_3)
% \d+ means "one or more digits"
pattern = '^\d+_\d+\.tif$'; 

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
    
    % Get image dimensions for the new border logic
    [imgHeight, imgWidth] = size(img_norm);
    
    % Segment
    cellLabel = segmentCells2D(cp, img_adj, ImageCellDiameter=image_cell_diameter);
    props = regionprops(cellLabel, img_norm, 'Area', 'MajorAxisLength', ...
        'MinorAxisLength', 'PixelIdxList', 'Centroid', 'Orientation', 'Circularity');
    
    thicknessMap = zeros(size(img_norm));
    I_bg = median(img_norm(cellLabel == 0));
    csvData = table();
    
    % Classification, Volume, and Table Logic
    for i = 1:length(props)
        aspectRatio = props(i).MajorAxisLength / props(i).MinorAxisLength;
        idx = props(i).PixelIdxList;
        
        % --- NEW LOGIC: 10-Pixel Centroid Border Exclusion ---
        % Centroid(1) is the X coordinate (width), Centroid(2) is the Y coordinate (height)
        cx = props(i).Centroid(1);
        cy = props(i).Centroid(2);
        
        isBorder = (cx <= 10) || (cx >= (imgWidth - 10)) || ...
                   (cy <= 10) || (cy >= (imgHeight - 10));
        
        currentVolume = 0; 
        
        if isBorder
            status = 'Border Exclusion';
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
    fig3D = figure('Name', ['3D Interactive - ', tif_files(k).name], 'Color', 'k'); 
    ax = axes('Parent', fig3D, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); 
    [X, Y] = meshgrid(1:size(img_norm,2), 1:size(img_norm,1));
    
    surf(ax, X, Y, thicknessMap, 'EdgeColor', 'none');
    colormap(turbo); 
    cb = colorbar;
    cb.Color = 'w'; 
    cb.Label.String = 'Thickness (\mum)';
    cb.Label.Color = 'w';
    
    view(45, 35); 
    axis tight; 
    grid on; 
    set(ax, 'GridColor', [0.3 0.3 0.3]); 
    
    camlight headlight; 
    lighting gouraud;
    
    title(['3D Reconstruction: ', tif_files(k).name], 'Color', 'w');
    zlabel('Height (\mum)', 'Color', 'w');
    pbaspect([1 1 0.5]);
    
    rotate3d on; 
end 
fprintf('Done. All images processed.\n');