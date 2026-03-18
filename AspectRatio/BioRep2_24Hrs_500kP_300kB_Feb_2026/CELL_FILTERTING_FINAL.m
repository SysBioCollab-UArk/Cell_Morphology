%% TIBD Morphological Analysis - Comprehensive Data Export
clear, clc, close all

% ------------------ Step 1: Configuration ------------------
main_folder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep2_24Hrs_500kP_300kB_Feb_2026';
outDir = fullfile(main_folder, 'Cluster_Results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cp = cellpose(model="cyto2");

% ------------------ Step 2: Nested Cluster Loop (Batch of 5x5) ------------------
for group = 1:1
    for imgNum = 1:5
        fileName = sprintf('%d_%d.tif', group, imgNum);
        fullPath = fullfile(main_folder, fileName);
        
        if ~exist(fullPath, 'file')
            fprintf('Skipping %s (file not found)\n', fileName);
            continue;
        end
        
        fprintf('Processing: %s\n', fileName);
        
        % Load and Pre-process
        img_raw = imread(fullPath);
        img = mean(single(img_raw), 3);
        img_adj = adapthisteq(img./65535);
        
        % --- SEGMENTATION & BORDER DETECTION ---
        allLabels = segmentCells2D(cp, img_adj, ImageCellDiameter=40);
        
        % Identify which labels touch the border
        clearedLabels = imclearborder(allLabels);
        borderMask = (allLabels > 0) & (clearedLabels == 0);
        borderIDs = unique(allLabels(borderMask)); 
        
        % Measure all requested parameters
        props = regionprops(allLabels, img, 'Area', 'MajorAxisLength', ...
            'MinorAxisLength', 'Orientation', 'Circularity', 'Centroid');
        
        if isempty(props), continue; end
        
        % Initialize Table with predefined Variable Names
        csvData = table();
        
        % ------------------ Step 3: Filtering & Classification ------------------
        for i = 1:length(props)
            % Calculate Aspect Ratio
            aspectRatio = props(i).MajorAxisLength / props(i).MinorAxisLength;
            
            % Check exclusion reasons
            isBorder = ismember(i, borderIDs);
            passMorphology = (aspectRatio > 1.1 && props(i).Circularity < 0.98);
            
            if isBorder
                status = 'Border Exclusion';
                textColor = 'red';
            elseif ~passMorphology
                status = 'Dead/Debris';
                textColor = 'red';
            else
                status = 'Viable Cells';
                textColor = 'green';
            end
            
            % Store color for the image labeling step
            props(i).DisplayColor = textColor;
            
            % Construct Data Row
            newRow = table(i, props(i).Area, props(i).MajorAxisLength, ...
                props(i).MinorAxisLength, aspectRatio, props(i).Orientation, ...
                props(i).Circularity, props(i).Centroid(1), props(i).Centroid(2), ...
                {status}, 'VariableNames', ...
                {'Cell_Number', 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
                'AspectRatio', 'Orientation', 'Circularity', 'Centroid_X', ...
                'Centroid_Y', 'Status'});
            
            csvData = [csvData; newRow];
        end
        
        % ------------------ Step 4: Save CSV & Labeled Image ------------------
        % 1. Save CSV for this image
        csvName = fullfile(outDir, [fileName(1:end-4), '_data.csv']);
        writetable(csvData, csvName);
        
        % 2. Create and Save Visual Check
        labelImg = uint8(img_adj .* 255); 
        labelImg = repmat(labelImg, 1, 1, 3); 
        
        for i = 1:length(props)
            if ~isempty(props(i).Centroid) && all(isfinite(props(i).Centroid))
                pos = props(i).Centroid;
                cellID = num2str(i);
                
                labelImg = insertText(labelImg, pos, cellID, ...
                    'AnchorPoint', 'Center', 'FontSize', 22, ...
                    'TextColor', props(i).DisplayColor, 'BoxOpacity', 0);
            end
        end
        
        imwrite(labelImg, fullfile(outDir, [fileName(1:end-4), '_labeledCells.tif']));
    end 
end 

fprintf('\nDone! CSVs and images generated for all batches.\n');