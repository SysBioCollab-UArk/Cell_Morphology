clear, clc, close all
% ------------------ Step 1: Set main folder ------------------
main_folder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep2_24Hrs_500kP_300kB_Feb_2026';

% Initialize Cellpose
cp = cellpose(model="cyto2");

% ------------------ Step 2: Control Loop ------------------
% Change these ranges to control which groups/images process
for group = 1:1
    for imgNum = 1:5
        
        fileName = sprintf('%d_%d.tif', group, imgNum);
        fullPath = fullfile(main_folder, fileName);
        
        if ~exist(fullPath, 'file')
            fprintf('Skipping %s (file not found)\n', fileName);
            continue;
        end
        
        fprintf('Processing: %s\n', fileName);
        img_raw = imread(fullPath);

        % Convert to grayscale if multi-channel
        if ndims(img_raw) == 3
            img = mean(single(img_raw), 3);   
        else
            img = single(img_raw);
        end

        % Normalize intensity (0 to 1 range for adapthisteq)
        img = img ./ 65535;         
        img_adj = adapthisteq(img);

        %% Segment cells 
        cellLabel = segmentCells2D(cp, img_adj, ImageCellDiameter=40);
        bw = cellLabel > 0;

        %% Get intensity 
        I_bg = median(img(cellLabel == 0));
        numCells = max(cellLabel(:));
        cellIntensity = zeros(1, numCells);
        for i = 1:numCells
            mask = (cellLabel == i);
            cellIntensity(i) = mean(img(mask));
        end
        I_cell = max(cellIntensity);

        %% Thickness reconstruction
        h_max = 12;   
        thickness = h_max .* (img - I_bg) / (I_cell - I_bg);
        thickness(thickness < 0) = 0;
        thickness(~bw) = 0; 

        %% Volume estimation
        cellVolume = zeros(1, numCells);
        for i = 1:numCells
            mask = (cellLabel == i);
            cellVolume(i) = sum(thickness(mask));
        end

        %% Visualization (Restored!)
        [X, Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
        figure('Color', 'k');
        surf(X, Y, thickness, 'EdgeColor', 'none');
        colormap turbo 
        colorbar;
        view(45, 30);
        axis tight; grid on;
        xlabel('X (px)'); ylabel('Y (px)'); zlabel('Thickness (µm)');
        title(['3D Map: ', fileName]);
        
        % Optional: Pause briefly to see the image before the next one pops up
        drawnow; 

        %% Extract quantitative thickness metrics
        mean_thickness = mean(thickness(bw));
        max_thickness  = max(thickness(bw));
        total_volume   = sum(thickness(bw)); 
        fprintf('--- Results for %s ---\n', fileName);
        fprintf('Mean thickness: %.2f µm | Max: %.2f µm\n', mean_thickness, max_thickness);

        %% Save the results as 32-bit Floating Point TIFF
        output_name = fullfile(main_folder, sprintf('thickness_%d_%d.tif', group, imgNum));
        t = Tiff(output_name, 'w');
        tagstruct.ImageLength     = size(thickness, 1);
        tagstruct.ImageWidth      = size(thickness, 2);
        tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample   = 32; 
        tagstruct.SamplesPerPixel = 1;
        tagstruct.RowsPerStrip    = 16;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SampleFormat    = Tiff.SampleFormat.IEEEFP; 
        tagstruct.Compression     = Tiff.Compression.None;
        t.setTag(tagstruct);
        t.write(single(thickness));
        t.close();
        
    end
end