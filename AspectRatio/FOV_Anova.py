# TIBD Automated Batch Processing - FOV and Sample-Level Analysis
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

cp = cellpose(model = "cyto2");

allCellData = table();

fprintf('Metadata loaded! Starting analysis for %d samples...\n', height(metaData));

% ------------------ Step 2: Cell-Level Processing ------------------

for i = 1:height(metaData)

    sampleID = metaData.Samples(i);
    currentCond = string(metaData.Condition{i});
    currentLine = string(metaData.CellLine{i});

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

% ------------------ Step 3: Save Cell-Level Data ------------------

if isempty(allCellData)
    error('No cells were detected by Cellpose.');
end

allCellData.Condition = categorical(string(allCellData.Condition));
allCellData.CellLine = categorical(string(allCellData.CellLine));

cellOutputFile = fullfile(mainFolder, 'BioRep1_CellLevel_Unfiltered_Results.csv');
writetable(allCellData, cellOutputFile);

fprintf('\nCell-level data saved to:\n%s\n', cellOutputFile);

% ------------------ Step 4: FOV-Level Summary ------------------

fovSummary = groupsummary( ...
    allCellData, ...
    {'Sample', 'Condition', 'CellLine', 'FOV'}, ...
    {'mean', 'median', 'std'}, ...
    'AspectRatio');

fovSummary.Properties.VariableNames{'mean_AspectRatio'} = 'Mean_AR';
fovSummary.Properties.VariableNames{'median_AspectRatio'} = 'Median_AR';
fovSummary.Properties.VariableNames{'std_AspectRatio'} = 'SD_AR';

fovOutputFile = fullfile(mainFolder, 'BioRep1_FOVLevel_Summary.csv');
writetable(fovSummary, fovOutputFile);

fprintf('FOV-level summary saved to:\n%s\n', fovOutputFile);

% ------------------ Step 5: Sample-Level Summary ------------------

sampleSummary = groupsummary( ...
    fovSummary, ...
    {'Sample', 'Condition', 'CellLine'}, ...
    {'mean', 'median', 'std'}, ...
    'Mean_AR');

sampleSummary.Properties.VariableNames{'mean_Mean_AR'} = 'Sample_Mean_AR';
sampleSummary.Properties.VariableNames{'median_Mean_AR'} = 'Sample_Median_AR';
sampleSummary.Properties.VariableNames{'std_Mean_AR'} = 'Sample_SD_AR';

sampleOutputFile = fullfile(mainFolder, 'BioRep1_SampleLevel_Summary.csv');
writetable(sampleSummary, sampleOutputFile);

fprintf('Sample-level summary saved to:\n%s\n', sampleOutputFile);

% ------------------ Step 6: Two-Way ANOVA Using FOV-Level Data ------------------

fprintf('\nRunning Two-Way ANOVA using FOV-level means...\n');

[p_fov, tbl_fov, stats_fov] = anovan( ...
    fovSummary.Mean_AR, ...
    {fovSummary.CellLine, fovSummary.Condition}, ...
    'model', 2, ...
    'varnames', {'CellLine', 'Condition'});

figure;
multcompare(stats_fov, 'Dimension', [1 2]);
title('Post-hoc Comparisons: FOV-Level Means');

% ------------------ Step 7: Two-Way ANOVA Using Sample-Level Data ------------------

fprintf('\nRunning Two-Way ANOVA using sample-level means...\n');

[p_sample, tbl_sample, stats_sample] = anovan( ...
    sampleSummary.Sample_Mean_AR, ...
    {sampleSummary.CellLine, sampleSummary.Condition}, ...
    'model', 2, ...
    'varnames', {'CellLine', 'Condition'});

figure;
multcompare(stats_sample, 'Dimension', [1 2]);
title('Post-hoc Comparisons: Sample-Level Means');

% ------------------ Step 8: Visualization - FOV-Level ------------------

figure('Color', 'w', 'Name', 'FOV-Level Cell Morphology Distribution');

boxchart(fovSummary.Condition, fovSummary.Mean_AR, ...
    'GroupByColor', fovSummary.CellLine);

grid on;
ylabel('Mean Aspect Ratio per FOV');
xlabel('Scaffold Condition');
title('FOV-Level Cell Elongation across BioRep1');
legend('Location', 'northeastoutside');
ylim([1, 4]);

% ------------------ Step 9: Visualization - Sample-Level ------------------

figure('Color', 'w', 'Name', 'Sample-Level Cell Morphology Distribution');

boxchart(sampleSummary.Condition, sampleSummary.Sample_Mean_AR, ...
    'GroupByColor', sampleSummary.CellLine);

grid on;
ylabel('Mean Aspect Ratio per Sample');
xlabel('Scaffold Condition');
title('Sample-Level Cell Elongation across BioRep1');
legend('Location', 'northeastoutside');
ylim([1, 4]);

fprintf('\nAnalysis complete.\n');