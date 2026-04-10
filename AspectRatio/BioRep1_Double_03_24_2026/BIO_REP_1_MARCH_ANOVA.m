%% TIBD Automated Batch Processing - Single Folder Edition
clear, clc, close all

% ------------------ Step 1: Configuration ------------------
% The single folder containing everything (Images + Excel)
mainFolder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026';

% Load the Sample ID file
metaData = readtable(fullfile(mainFolder, 'SampleID.xlsx'));

% Initialize Cellpose and Results Table
cp = cellpose(model="cyto2");
allAnalysisData = table(); 

fprintf('Metadata loaded! Starting analysis for %d samples...\n', height(metaData));

% ------------------ Step 2: Processing Loop ------------------
for i = 1:height(metaData)
    sampleID = metaData.Samples(i);      % Column "Samples"
    currentCond = metaData.Condition{i}; % Column "Condition"
    currentLine = metaData.CellLine{i};  % Column "CellLine"
    
    for m = 1:5
        fileName = sprintf('%d_%d.tif', sampleID, m);
        fullPath = fullfile(mainFolder, fileName);
        
        if ~exist(fullPath, 'file')
            continue; % Skip if the image isn't there
        end
        
        fprintf('Processing: %s [%s | %s]\n', fileName, currentCond, currentLine);
        
        % --- Image Analysis ---
        img = mean(single(imread(fullPath)), 3);
        img_adj = adapthisteq(img./65535);
        allLabels = segmentCells2D(cp, img_adj, ImageCellDiameter=40);
        
        % Border Exclusion
        cleared = imclearborder(allLabels);
        borderIDs = unique(allLabels(allLabels > 0 & cleared == 0));
        
        props = regionprops(allLabels, img, 'MajorAxisLength', 'MinorAxisLength', 'Circularity');
        if isempty(props), continue; end
        
        % --- Apply Your Exclusion Logic ---
        tempTable = table();
        for k = 1:length(props)
            aspectRatio = props(k).MajorAxisLength / props(k).MinorAxisLength;
            isBorder = ismember(k, borderIDs);
            passMorphology = (aspectRatio > 1.1 && props(k).Circularity < 0.98);
            
            if ~isBorder && passMorphology
                newRow = table(sampleID, {currentCond}, {currentLine}, aspectRatio, ...
                    'VariableNames', {'Sample', 'Condition', 'CellLine', 'AspectRatio'});
                tempTable = [tempTable; newRow];
            end
        end
        allAnalysisData = [allAnalysisData; tempTable];
    end
end

% ------------------ Step 3: Statistics ------------------
if ~isempty(allAnalysisData)
    % Format for ANOVA
    allAnalysisData.Condition = categorical(string(allAnalysisData.Condition));
    allAnalysisData.CellLine = categorical(string(allAnalysisData.CellLine));
    
    % Two-Way ANOVA
    [p, tbl, stats] = anovan(allAnalysisData.AspectRatio, ...
        {allAnalysisData.CellLine, allAnalysisData.Condition}, ...
        'model', 2, 'varnames', {'CellLine', 'Condition'});
    
    % Visual Comparison Plot
    figure;
    multcompare(stats, 'Dimension', [1 2]);
    
    % Save results
    writetable(allAnalysisData, fullfile(mainFolder, 'BioRep1_Final_Results.csv'));
    fprintf('\nSuccess! Results saved to BioRep1_Final_Results.csv\n');
else
    disp('Check failed: No viable cells were found in the images.');
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