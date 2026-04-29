%% TIBD Unified High-Precision Analysis - BioRep2 (v3 - Phase Contrast Aware)
clear, clc, close all

% ------------------ Step 1: Configuration ------------------
%mainFolder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/2026_04_09_10_38_48--BioRep2_Double_Cells_04_09_2026_Final';
mainFolder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026';
%metaData   = readtable(fullfile(mainFolder, 'ID_Samples_BioRep2_Double.xlsx'));
metaData   = readtable(fullfile(mainFolder, 'SampleID.xlsx'));
outDir = fullfile(mainFolder, 'High_Precision_Analysis');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('Initializing Cellpose Model...\n');
cp = cellpose(model="cyto2");

allAnalysisData = table();
fprintf('Metadata loaded. Starting analysis for %d samples...\n', height(metaData));

% ============================================================
% HELPER: compute per-cell phase-contrast-aware features
%   img      = raw single-precision image (0–65535 range)
%   allLabels = integer label matrix from Cellpose
%   Returns a struct array with one entry per label index
% ============================================================
function feat = computePhaseFeatures(img, allLabels, nLabels)
    feat(nLabels).intensityStd  = 0;
    feat(nLabels).haloRatio     = 1;
    feat(nLabels).maxIntensity  = 0;
    feat(nLabels).perimRough    = 1;
    feat(nLabels).localContrast = 0;

    se_erode  = strel('disk', 4);   % erosion disk for interior
    se_dilate = strel('disk', 8);   % dilation disk for local background

    for k = 1:nLabels
        mask = (allLabels == k);
        if ~any(mask(:)), continue; end
        pixels = img(mask);

        % ------ Feature 1: Intensity std dev -----------------
        feat(k).intensityStd = std(single(pixels));

        % ------ Feature 2: Halo ratio ------------------------
        interior = imerode(mask, se_erode);
        boundary = mask & ~interior;
        if sum(interior(:)) > 10 && sum(boundary(:)) > 5
            feat(k).haloRatio = mean(img(boundary)) / (mean(img(interior)) + 1);
        else
            feat(k).haloRatio = 1;
        end

        % ------ Feature 3: Max pixel intensity ---------------
        feat(k).maxIntensity = max(pixels);

        % ------ Feature 4: Perimeter roughness ---------------
        perim = sum(bwperim(mask, 8), 'all');
        area  = sum(mask, 'all');
        if area > 0
            feat(k).perimRough = (double(perim)^2) / (4 * pi * double(area));
        else
            feat(k).perimRough = 1;
        end

        % ------ Feature 5: Local contrast -------------------
        localRegion = imdilate(mask, se_dilate) & ~mask;
        if sum(localRegion(:)) > 20
            feat(k).localContrast = mean(pixels) - mean(img(localRegion));
        else
            feat(k).localContrast = 0;
        end
    end
end

% ------------------ Step 2: Processing Loop ------------------
for i = 1:height(metaData)
    sampleID     = metaData.Samples(i);
    currentCond  = metaData.Condition{i};
    currentLine  = metaData.CellLine{i};

    for m = 1:5
        fileName = sprintf('%d_%d.tif', sampleID, m);
        fullPath = fullfile(mainFolder, fileName);
        if ~exist(fullPath, 'file'), continue; end

        fprintf('Processing: %s [%s | %s]\n', fileName, currentCond, currentLine);

        % ---- Pre-processing ----
        img_raw = imread(fullPath);
        img     = mean(single(img_raw), 3);           % raw 16-bit for measurements
        img_adj = adapthisteq(img ./ 65535, ...
                      'ClipLimit', 0.015);            % CLAHE only for segmentation
                      
        % NEW: Get image dimensions for the 10-pixel border rule
        [imgHeight, imgWidth] = size(img);

        % ---- Segmentation ----
        allLabels = segmentCells2D(cp, img_adj, ImageCellDiameter=40);

        % ---- Standard regionprops ----
        props = regionprops(allLabels, img, ...
            'Area', 'MajorAxisLength', 'MinorAxisLength', ...
            'Orientation', 'Circularity', 'Centroid', ...
            'Solidity', 'MeanIntensity', 'Eccentricity');

        if isempty(props), continue; end
        nLabels = length(props);

        % ================================================================
        % COMPUTE PHASE-CONTRAST-AWARE FEATURES (the main new addition)
        % ================================================================
        feat = computePhaseFeatures(img, allLabels, nLabels);

        % ----------------------------------------------------------------
        % POPULATION-LEVEL THRESHOLDS
        % ----------------------------------------------------------------
        allMeans  = [props.MeanIntensity];
        allStds   = [feat.intensityStd];
        allHalos  = [feat.haloRatio];
        allMaxes  = [feat.maxIntensity];

        meanCutoff    = median(allMeans) + 2.5 * mad(allMeans, 1);
        hardMeanCut   = median(allMeans) + 5.0 * mad(allMeans, 1);
        stdCutoff     = median(allStds)  + 2.0 * mad(allStds,  1);
        haloCutoff    = median(allHalos) + 2.5 * mad(allHalos, 1);
        maxCutoff     = prctile(allMaxes, 95);

        tempTable = table();

        for k = 1:nLabels
            if props(k).MinorAxisLength > 0
                aspectRatio = props(k).MajorAxisLength / props(k).MinorAxisLength;
            else
                aspectRatio = 1;
            end
            
            % NEW 10-Pixel Centroid Rule
            % Centroid(1) is X (width), Centroid(2) is Y (height)
            cx = props(k).Centroid(1);
            cy = props(k).Centroid(2);
            is_border = (cx <= 10) || (cx >= imgWidth - 10) || ...
                        (cy <= 10) || (cy >= imgHeight - 10);

            % ---- Compute boolean flags ----
            is_debris           = props(k).Area < 150;
            
            % MEAN-BASED
            is_too_bright       = props(k).MeanIntensity > meanCutoff;
            is_extremely_bright = props(k).MeanIntensity > hardMeanCut;
            
            is_high_variance    = feat(k).intensityStd > stdCutoff;
            is_ghost_halo       = feat(k).haloRatio > haloCutoff;
            is_max_outlier      = feat(k).maxIntensity > maxCutoff;
            is_blebbing         = feat(k).perimRough > 1.35 && is_too_bright;
            is_locally_bright   = feat(k).localContrast > (0.15 * 65535);

            % Shape flags 
            is_round_dead       = props(k).Circularity > 0.72;
            is_viable_size      = props(k).Area >= 200;

            % ================================================================
            % DECISION TREE
            % ================================================================
            if is_border
                status     = 'Excluded_Border';
                textColor  = [1 1 0];           % yellow
            elseif is_debris
                status     = 'Excluded_Debris';
                textColor  = [1 0 0];
            elseif is_extremely_bright
                status     = 'Excluded_Dead_Extreme';
                textColor  = [0.7 0 0];
            elseif is_ghost_halo && is_round_dead
                status     = 'Excluded_Dead_Ghost';
                textColor  = [1 0.2 0];
            elseif is_high_variance && is_round_dead
                status     = 'Excluded_Dead_Ghost';
                textColor  = [1 0.2 0];
            elseif is_too_bright && is_round_dead
                status     = 'Excluded_Dead_Bright';
                textColor  = [1 0 0];
            elseif is_blebbing
                status     = 'Excluded_Dead_Blebbing';
                textColor  = [1 0.4 0];
            elseif is_max_outlier && is_locally_bright
                status     = 'Excluded_Dead_Fragment';
                textColor  = [0.9 0.2 0.2];
            elseif is_too_bright && props(k).Area < 500
                status     = 'Excluded_Debris_Bright';
                textColor  = [1 0 0];
            elseif is_viable_size
                % *** REVISED VIABLE GATE ***
                % If it is large enough and did NOT trigger any of the dead, ghost, 
                % blebbing, or border flags above, we keep it! 
                % No pre-filtering for shape to prevent confirmation bias.
                status     = 'Viable';
                textColor  = [0 1 0];
            else
                % Borderline — flag orange for manual QC review
                status     = 'Excluded_Other';
                textColor  = [1 0.6 0];
            end

            newRow = table(sampleID, {currentCond}, {currentLine}, aspectRatio, ...
                props(k).Orientation, props(k).Area, props(k).Circularity, ...
                props(k).Solidity, props(k).MeanIntensity, props(k).Eccentricity, ...
                feat(k).intensityStd, feat(k).haloRatio, feat(k).maxIntensity, ...
                feat(k).perimRough, feat(k).localContrast, {status}, ...
                'VariableNames', {'Sample','Condition','CellLine','AspectRatio', ...
                'Orientation','Area','Circularity','Solidity','Intensity','Eccentricity', ...
                'IntensityStd','HaloRatio','MaxIntensity','PerimRoughness', ...
                'LocalContrast','Status'});
                
            tempTable = [tempTable; newRow];
            props(k).DisplayColor = textColor;
        end
        
        allAnalysisData = [allAnalysisData; tempTable];

        % Visual overlay
        labelImg = repmat(uint8(img_adj .* 255), 1, 1, 3);
        for k = 1:length(props)
            if isfinite(props(k).Centroid(1))
                labelImg = insertText(labelImg, props(k).Centroid, num2str(k), ...
                    'AnchorPoint','Center','FontSize',14, ...
                    'TextColor', props(k).DisplayColor, 'BoxOpacity', 0);
            end
        end
        imwrite(labelImg, fullfile(outDir, [fileName(1:end-4), '_Processed.tif']));
    end
end

% ------------------ Step 3: Statistics & ANOVA ------------------
if ~isempty(allAnalysisData)
    writetable(allAnalysisData, fullfile(mainFolder, 'BioRep2_HighPrecision_ALL_DATA.csv'));
    
    viableData = allAnalysisData(strcmp(allAnalysisData.Status, 'Viable'), :);
    
    if ~isempty(viableData)
        viableData.Condition = removecats(categorical(string(viableData.Condition)));
        viableData.CellLine  = removecats(categorical(string(viableData.CellLine)));
        
        [~, ~, stats] = anovan(viableData.AspectRatio, ...
            {viableData.CellLine, viableData.Condition}, ...
            'model', 2, 'varnames', {'CellLine','Condition'});
            
        figure; multcompare(stats, 'Dimension', [1 2]);
        title('Multcompare: Viable Cell Aspect Ratio');
        
        figure('Color','w','Name','Cell Elongation Distribution');
        boxchart(viableData.Condition, viableData.AspectRatio, ...
                 'GroupByColor', viableData.CellLine);
        ylabel('Aspect Ratio (L/W)'); xlabel('Condition');
        title('Viable Cell Morphology (Filtered)'); grid on; legend;
    end

    % ---- Classification audit printout ----
    fprintf('\n--- Cell Classification Summary ---\n');
    allStatuses = string(allAnalysisData.Status);
    cats = unique(allStatuses);
    for s = 1:length(cats)
        n = sum(allStatuses == cats(s));
        pct = 100 * n / height(allAnalysisData);
        fprintf('  %-35s : %4d  (%.1f%%)\n', cats(s), n, pct);
    end
    fprintf('  %-35s : %4d\n', 'TOTAL', height(allAnalysisData));
else
    disp('No objects detected.');
end

fprintf('\nAnalysis complete.\n');