%% TIBD Automated Batch Processing - Cell-Level Six-Group Analysis
clear; clc; close all;

% ------------------ Step 1: Configuration ------------------

mainFolder = '/Users/alexandravega/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026';
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

    end
end

% ------------------ Step 4: Save Cell-Level Data ------------------

if isempty(allCellData)
    error('No cells were detected by Cellpose.');
end

allCellData.Condition = categorical(strtrim(string(allCellData.Condition)));
allCellData.CellLine = categorical(strtrim(string(allCellData.CellLine)));

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

outputFile = fullfile(mainFolder, 'BioRep1_CellLevel_SixGroup_Results.csv');
writetable(allCellData, outputFile);

fprintf('\nCell-level data saved to:\n%s\n', outputFile);

% ------------------ Step 5: Box Plot ------------------

% figure('Color', 'w', 'Name', 'Boxplot - Aspect Ratio');
% 
% boxchart(allCellData.Group, allCellData.AspectRatio);
% 
% grid on;
% ylabel('Aspect Ratio');
% xlabel('Experimental Group');
% title('Cell Aspect Ratio Across Six Experimental Conditions');
% ylim([1, 4]);


% ------------------ Step 5: Box Plot (No Outliers) ------------------

figure('Color', 'w', 'Name', 'Boxplot - Aspect Ratio');

b = boxchart(allCellData.Group, allCellData.AspectRatio);

% Remove outlier markers (visual only)
b.MarkerStyle = 'none';

grid on;
ylabel('Aspect Ratio');
xlabel('Experimental Group');
title('Cell Aspect Ratio Across Six Experimental Conditions');

ylim([0, 4]);   % Lower limit = 0

% ------------------ Step 6: Violin Plot ------------------
% Requires violinplot.m to be available in your MATLAB path.

figure('Color', 'w', 'Name', 'Violin Plot - Aspect Ratio');

violinplot(allCellData.Group, allCellData.AspectRatio);

grid on;
ylabel('Aspect Ratio');
xlabel('Experimental Group');
title('Cell Aspect Ratio Distribution Across Six Experimental Conditions');
ylim([1, 4]);

% ------------------ Step 7: Targeted Statistical Comparisons ------------------

% Remove cells with undefined aspect ratios
validData = allCellData(~isnan(allCellData.AspectRatio), :);

Bone_Control = validData.AspectRatio(validData.Group == "Bone Clone-TC");
Bone_Random = validData.AspectRatio(validData.Group == "Bone Clone-Random");
Bone_Aligned = validData.AspectRatio(validData.Group == "Bone Clone-Aligned");

Parental_Control = validData.AspectRatio(validData.Group == "Parental-TC");
Parental_Random = validData.AspectRatio(validData.Group == "Parental-Random");
Parental_Aligned = validData.AspectRatio(validData.Group == "Parental-Aligned");

fprintf('\n===== SAMPLE SIZE PER GROUP =====\n');
fprintf('Bone Clone-TC: %d cells\n', numel(Bone_Control));
fprintf('Bone Clone-Random: %d cells\n', numel(Bone_Random));
fprintf('Bone Clone-Aligned: %d cells\n', numel(Bone_Aligned));
fprintf('Parental-TC: %d cells\n', numel(Parental_Control));
fprintf('Parental-Random: %d cells\n', numel(Parental_Random));
fprintf('Parental-Aligned: %d cells\n', numel(Parental_Aligned));

fprintf('\n===== Mann-Whitney U Tests / ranksum =====\n');

p_Bone_CR = runRanksumSafe(Bone_Control, Bone_Random, ...
    'Bone Clone: TC vs Random');

p_Bone_CA = runRanksumSafe(Bone_Control, Bone_Aligned, ...
    'Bone Clone: TC vs Aligned');

p_Bone_RA = runRanksumSafe(Bone_Random, Bone_Aligned, ...
    'Bone Clone: Random vs Aligned');

p_Parental_CR = runRanksumSafe(Parental_Control, Parental_Random, ...
    'Parental: TC vs Random');

p_Parental_CA = runRanksumSafe(Parental_Control, Parental_Aligned, ...
    'Parental: TC vs Aligned');

p_Parental_RA = runRanksumSafe(Parental_Random, Parental_Aligned, ...
    'Parental: Random vs Aligned');

% ------------------ Step 8: Save Statistical Results ------------------

comparisonNames = { ...
    'Bone Clone TC vs Random'; ...
    'Bone Clone TC vs Aligned'; ...
    'Bone Clone Random vs Aligned'; ...
    'Parental TC vs Random'; ...
    'Parental TC vs Aligned'; ...
    'Parental Random vs Aligned' ...
};

pValues = [ ...
    p_Bone_CR; ...
    p_Bone_CA; ...
    p_Bone_RA; ...
    p_Parental_CR; ...
    p_Parental_CA; ...
    p_Parental_RA ...
];

% Bonferroni correction
pValues_adj = min(pValues * 6, 1);

statsTable = table( ...
    comparisonNames, ...
    pValues, ...
    pValues_adj, ...
    'VariableNames', {'Comparison', 'P_Value', 'Adjusted_P'} ...
);

statsOutputFile = fullfile(mainFolder, 'BioRep1_Targeted_Stats_CellLevel.csv');
writetable(statsTable, statsOutputFile);

fprintf('\nTargeted statistical results saved to:\n%s\n', statsOutputFile);

fprintf('\nAnalysis complete.\n');

% ------------------ Local Function ------------------
% This protects the script from crashing if a group is empty.

function pValue = runRanksumSafe(group1, group2, comparisonName)

    group1 = group1(~isnan(group1));
    group2 = group2(~isnan(group2));

    if isempty(group1) || isempty(group2)
        fprintf('%s skipped: one or both groups are empty.\n', comparisonName);
        pValue = NaN;
    else
        [pValue, ~] = ranksum(group1, group2);
        fprintf('%s p = %.5f\n', comparisonName, pValue);
    end

end