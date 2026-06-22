%% TIBD Automated Batch Processing - Cell-Level Six Group Analysis
clear; clc; close all;

% ------------------ Step 1: Configuration ------------------

mainFolder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026';
metadataFile = fullfile(mainFolder, 'SampleID.xlsx');

if ~exist(metadataFile, 'file')
    error('Metadata file not found: %s', metadataFile);
end

metaData = readtable(metadataFile);

requiredColumns = {'Samples', 'Condition', 'CellLine'};
if ~all(ismember(requiredColumns, metaData.Properties.VariableNames))
    error('Metadata file must contain columns: Samples, Condition, CellLine');
end

% ------------------ Step 2: Initialize Cellpose ------------------

cp = cellpose(model = "cyto2");
allCellData = table();

fprintf('Metadata loaded! Starting cell-level analysis for %d samples...\n', height(metaData));

% ------------------ Step 3: Process Images ------------------

for i = 1:height(metaData)

    sampleID = metaData.Samples(i);
    currentCond = strtrim(string(metaData.Condition{i}));
    currentLine = strtrim(string(metaData.CellLine{i}));

    for m = 1:5

        fileName = sprintf('%d_%d.tif', sampleID, m);
        fullPath = fullfile(mainFolder, fileName);

        if ~exist(fullPath, 'file')
            fprintf('File not found, skipping: %s\n', fileName);
            continue;
        end

        fprintf('Processing: %s [%s | %s]\n', fileName, currentCond, currentLine);

        imgRaw = imread(fullPath);
        img = mean(single(imgRaw), 3);

        imgNorm = img ./ max(img(:));
        img_adj = adapthisteq(imgNorm);

        allLabels = segmentCells2D(cp, img_adj, ImageCellDiameter = 40);

        props = regionprops(allLabels, img, ...
            'Area', 'MajorAxisLength', 'MinorAxisLength', 'Circularity');

        if isempty(props)
            fprintf('No cells detected in: %s\n', fileName);
            continue;
        end

        tempTable = table();

        for k = 1:length(props)

            if props(k).MinorAxisLength == 0
                aspectRatio = NaN;
            else
                aspectRatio = props(k).MajorAxisLength / props(k).MinorAxisLength;
            end

            newRow = table( ...
                sampleID, ...
                {currentCond}, ...
                {currentLine}, ...
                m, ...
                k, ...
                aspectRatio, ...
                props(k).Area, ...
                props(k).Circularity, ...
                'VariableNames', {'Sample', 'Condition', 'CellLine', ...
                                  'FOV', 'CellID', 'AspectRatio', ...
                                  'Area', 'Circularity'} ...
            );

            tempTable = [tempTable; newRow];

        end

        allCellData = [allCellData; tempTable];

        % ------------------ Visualization: Overlay on Image ------------------
        % Every cell identified by Cellpose is outlined in green.
        % Boundaries are thickened using imdilate for solid, visible outlines.

        % Convert grayscale to uint8 RGB base image
        imgDisp = uint8(imgNorm * 255);
        imgRGB  = cat(3, imgDisp, imgDisp, imgDisp);

        % Build a single binary boundary mask for ALL labeled cells,
        % then dilate it to make lines thick and solid (like the reference image).
        boundaryMask = false(size(allLabels));
        allIDs = unique(allLabels(:));
        allIDs = allIDs(allIDs > 0);   % remove background

        for k = allIDs'
            cellMask = (allLabels == k);
            % erode the cell mask and XOR to get the boundary ring
            erodedMask = imerode(cellMask, strel('disk', 1));
            boundaryMask = boundaryMask | (cellMask & ~erodedMask);
        end

        % Dilate boundary mask to get thick, solid lines (adjust radius for line width)
        se = strel('disk', 2);   % radius 2 → ~5px wide line; increase for thicker
        thickBoundary = imdilate(boundaryMask, se);

        % Paint thick green boundary onto the RGB image
        imgRGB(:,:,1) = imgRGB(:,:,1) .* uint8(~thickBoundary);   % zero out R
        imgRGB(:,:,2) = imgRGB(:,:,2) .* uint8(~thickBoundary) + uint8(thickBoundary) * 220;  % set G
        imgRGB(:,:,3) = imgRGB(:,:,3) .* uint8(~thickBoundary);   % zero out B

        % Save overlay figure — title in black for readability
        hFig = figure('Color', 'w', 'Visible', 'off');
        imshow(imgRGB);
        t = title(sprintf('Green: Viable | Red: Excluded - %d_%d.tif', sampleID, m), ...
              'FontWeight', 'bold', 'FontName', 'Monospaced');
        t.Color = [0 0 0];   % black text

        overlayName = sprintf('%d_%d_overlay.png', sampleID, m);
        overlayPath = fullfile(mainFolder, overlayName);
        exportgraphics(hFig, overlayPath, 'Resolution', 150);
        close(hFig);

        fprintf('  Overlay saved: %s\n', overlayName);

    end
end

% ------------------ Step 4: Save Cell-Level Data ------------------

if isempty(allCellData)
    error('No cells were detected by Cellpose.');
end

allCellData.Condition = categorical(strtrim(string(allCellData.Condition)));
allCellData.CellLine  = categorical(strtrim(string(allCellData.CellLine)));

allCellData.Group = strcat(string(allCellData.CellLine), "-", string(allCellData.Condition));
allCellData.Group = categorical(allCellData.Group);

disp('Groups found in the data:');
disp(categories(allCellData.Group));

groupOrder = [ ...
    "Bone Clone-TC", ...
    "Bone Clone-Random", ...
    "Bone Clone-Aligned", ...
    "Parental-TC", ...
    "Parental-Random", ...
    "Parental-Aligned" ...
];

existingGroups = string(categories(allCellData.Group));
validGroupOrder = groupOrder(ismember(groupOrder, existingGroups));

missingGroups = groupOrder(~ismember(groupOrder, existingGroups));

if ~isempty(missingGroups)
    disp('Warning: These expected groups were NOT found in your data:');
    disp(missingGroups');
end

allCellData.Group = reordercats(allCellData.Group, validGroupOrder);

% Save to Excel (primary output)
excelFileName = fullfile(mainFolder, 'BioRep1_CellData.xlsx');
writetable(allCellData, excelFileName);
fprintf('Success! Data saved to Excel for Python analysis.\n');

fprintf('\nAnalysis complete.\n');
fprintf('  Excel file : %s\n', excelFileName);
fprintf('  Overlays   : saved alongside each source .tif in:\n             %s\n', mainFolder);