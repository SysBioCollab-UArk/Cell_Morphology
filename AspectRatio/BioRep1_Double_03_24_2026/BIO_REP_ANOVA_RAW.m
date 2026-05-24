%% TIBD Automated Batch Processing - NO EXCLUSIONS
clear, clc, close all

% ------------------ Step 1: Configuration ------------------
mainFolder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/2026_04_09_10_38_48--BioRep2_Double_Cells_04_09_2026_Final';
%mainFolder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026';

% Added the .xlsx extension (assumed) to ensure the table loads correctly
metaData = readtable(fullfile(mainFolder, 'ID_Samples_BioRep2_Double.xlsx'));
%metaData = readtable(fullfile(mainFolder, 'SampleID.xlsx'));   %BioRep 1

cp = cellpose(model="cyto2");
allAnalysisData = table(); 

fprintf('Metadata loaded! Starting RAW analysis for %d samples...\n', height(metaData));

% ------------------ Step 2: Processing Loop ------------------
for i = 1:height(metaData)
    sampleID = metaData.Samples(i);      
    currentCond = metaData.Condition{i}; 
    currentLine = metaData.CellLine{i};  
    
    for m = 1:5
        fileName = sprintf('%d_%d.tif', sampleID, m);
        fullPath = fullfile(mainFolder, fileName);
        
        if ~exist(fullPath, 'file'), continue; end
        
        fprintf('Processing Raw Data: %s [%s | %s]\n', fileName, currentCond, currentLine);
        
        % Load image
        img = mean(single(imread(fullPath)), 3);
        img_adj = adapthisteq(img./65535);
        
        % Run Segmentation
        allLabels = segmentCells2D(cp, img_adj, ImageCellDiameter=40);
        
        % --- UPDATE 1: Measure ALL parameters ---
        props = regionprops(allLabels, img, 'Area', 'MajorAxisLength', ...
            'MinorAxisLength', 'Circularity', 'Orientation', 'Solidity', ...
            'Eccentricity', 'Perimeter', 'Extent', 'MeanIntensity', 'MaxIntensity');
            
        if isempty(props), continue; end
        
        % Calculate global background intensity for Local Contrast math
        bg_mask = (allLabels == 0);
        I_bg = mean(img(bg_mask), 'omitnan');
        
        tempTable = table();
        for k = 1:length(props)
            % Native Geometric Metrics
            aspectRatio = props(k).MajorAxisLength / props(k).MinorAxisLength;
            
            % Native Intensity Metrics
            intensity = props(k).MeanIntensity;
            maxIntensity = props(k).MaxIntensity;
            
            % --- UPDATE 2: Custom Derived Metrics ---
            % Local Contrast (Cell intensity minus empty background intensity)
            localContrast = intensity - I_bg; 
            
            % Halo Ratio (How bright the edge halo is compared to main body)
            haloRatio = maxIntensity / intensity;
            
            % Perimeter Roughness (Deviation from a perfectly smooth circular edge)
            perimRoughness = (props(k).Perimeter^2) / (4 * pi * props(k).Area);
            
            % --- EXCLUSION CRITERIA REMOVED ---
            % --- UPDATE 3: Save EVERYTHING to the row ---
            newRow = table(sampleID, {currentCond}, {currentLine}, aspectRatio, ...
                props(k).Area, props(k).Circularity, props(k).Orientation, ...
                props(k).Solidity, props(k).Eccentricity, props(k).Perimeter, ...
                props(k).Extent, intensity, maxIntensity, localContrast, ...
                haloRatio, perimRoughness, ...
                'VariableNames', {'Sample', 'Condition', 'CellLine', 'AspectRatio', ...
                'Area', 'Circularity', 'Orientation', 'Solidity', ...
                'Eccentricity', 'Perimeter', 'Extent', 'Intensity', 'MaxIntensity', ...
                'LocalContrast', 'HaloRatio', 'PerimRoughness'});
            
            tempTable = [tempTable; newRow];
        end
        allAnalysisData = [allAnalysisData; tempTable];
    end
end

% ------------------ Step 3: Statistics & Saving ------------------
if ~isempty(allAnalysisData)
    % Prepare for ANOVA
    allAnalysisData.Condition = categorical(string(allAnalysisData.Condition));
    allAnalysisData.CellLine = categorical(string(allAnalysisData.CellLine));
    
    % Save full dataset before any stats
    %writetable(allAnalysisData, fullfile(mainFolder, 'BioRep2_Full_Unfiltered_Results.csv'));
    writetable(allAnalysisData, fullfile(mainFolder, 'BioRep1_Full_Unfiltered_Results.csv'));
    
    % Two-Way ANOVA on the raw data (Using AspectRatio as the default test case)
    [p, tbl, stats] = anovan(allAnalysisData.AspectRatio, ...
        {allAnalysisData.CellLine, allAnalysisData.Condition}, ...
        'model', 2, 'varnames', {'CellLine', 'Condition'});
    
    % Visualization
    figure;
    multcompare(stats, 'Dimension', [1 2]);
    
    fprintf('\nSuccess! All raw data saved to BioRep1_Full_Unfiltered_Results.csv\n');
else
    disp('Error: No cells were detected by Cellpose.');
end

% ------------------ Step 4: Advanced Visualization ------------------
figure('Color', 'w', 'Name', 'Cell Morphology Distribution');

% Create a boxchart grouped by Condition and colored by CellLine
bn = boxchart(allAnalysisData.Condition, allAnalysisData.AspectRatio, ...
    'GroupByColor', allAnalysisData.CellLine);

% Make it look professional
bn(1).SeriesIndex = 1; % Color for Bone Clone
bn(2).SeriesIndex = 2; % Color for Parental
grid on
ylabel('Aspect Ratio (Length/Width)');
xlabel('Scaffold Condition');
title('Distribution of Cell Elongation across BioRep1');
legend('Location', 'northeastoutside');

% Optional: Set Y-axis limits to focus on the bulk of the data
ylim([1, 4]);