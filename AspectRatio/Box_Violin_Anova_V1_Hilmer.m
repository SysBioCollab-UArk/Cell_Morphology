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

% ------------------ Step 5: Box Plot (No Outliers) ------------------
figure('Color', 'w', 'Name', 'Boxplot - Aspect Ratio');
b = boxchart(allCellData.Group, allCellData.AspectRatio);
b.MarkerStyle = 'none'; 
grid on;
ylabel('Aspect Ratio');
xlabel('Experimental Group');
% Title removed per request
ylim([0, 4]); % Lower limit set to 0

% ------------------ Step 6: Violin Plot ------------------
figure('Color', 'w', 'Name', 'Violin Plot - Aspect Ratio');
violinplot(allCellData.Group, allCellData.AspectRatio);
grid on;
ylabel('Aspect Ratio');
xlabel('Experimental Group');
% Title removed per request
ylim([0, 4]); % Lower limit set to 0

% ------------------ Step 7: Comprehensive Statistical Comparisons ------------------

validData = allCellData(~isnan(allCellData.AspectRatio), :);

Bone_Control     = validData.AspectRatio(validData.Group == "Bone Clone-TC");
Bone_Random      = validData.AspectRatio(validData.Group == "Bone Clone-Random");
Bone_Aligned     = validData.AspectRatio(validData.Group == "Bone Clone-Aligned");
Parental_Control = validData.AspectRatio(validData.Group == "Parental-TC");
Parental_Random  = validData.AspectRatio(validData.Group == "Parental-Random");
Parental_Aligned = validData.AspectRatio(validData.Group == "Parental-Aligned");

fprintf('\n===== SAMPLE SIZE PER GROUP =====\n');
fprintf('  Bone Clone-TC:      %d cells\n', numel(Bone_Control));
fprintf('  Bone Clone-Random:  %d cells\n', numel(Bone_Random));
fprintf('  Bone Clone-Aligned: %d cells\n', numel(Bone_Aligned));
fprintf('  Parental-TC:        %d cells\n', numel(Parental_Control));
fprintf('  Parental-Random:    %d cells\n', numel(Parental_Random));
fprintf('  Parental-Aligned:   %d cells\n', numel(Parental_Aligned));

% ---- 7a. Normality Testing (Lilliefors on log-transformed data per group) ----
fprintf('\n===== NORMALITY — Lilliefors on log(AspectRatio) =====\n');
fprintf('  (p < 0.05 → reject normality; if all "normal" then Welch t is trustworthy)\n');

groupLabels_norm = { ...
    'Bone Clone-TC',    'Bone Clone-Random', 'Bone Clone-Aligned', ...
    'Parental-TC',      'Parental-Random',   'Parental-Aligned'};
groupArrays_norm = { ...
    Bone_Control, Bone_Random, Bone_Aligned, ...
    Parental_Control, Parental_Random, Parental_Aligned};

normalityResults = table();
for g = 1:6
    gData   = groupArrays_norm{g};
    gData   = gData(~isnan(gData) & gData > 0);
    logData = log(gData);
    [h_lil, p_lil] = lillietest(logData);
    if h_lil
        normStr = 'NOT normal';
    else
        normStr = 'normal';
    end
    fprintf('  %-25s  p = %.5f  (%s)\n', groupLabels_norm{g}, p_lil, normStr);
    normalityResults = [normalityResults; table( ...
        {groupLabels_norm{g}}, p_lil, logical(h_lil), ...
        'VariableNames', {'Group', 'Lilliefors_P', 'RejectNormality'})];
end

% ---- 7b. Pairwise comparisons — all four test families ----
comparisons = { ...
    Bone_Control,     Bone_Random,     'Bone Clone: TC vs Random'; ...
    Bone_Control,     Bone_Aligned,    'Bone Clone: TC vs Aligned'; ...
    Bone_Random,      Bone_Aligned,    'Bone Clone: Random vs Aligned'; ...
    Parental_Control, Parental_Random,  'Parental: TC vs Random'; ...
    Parental_Control, Parental_Aligned, 'Parental: TC vs Aligned'; ...
    Parental_Random,  Parental_Aligned, 'Parental: Random vs Aligned'; ...
};

nComp    = size(comparisons, 1);
allStats = table();

fprintf('\n===== PAIRWISE TESTS — All Families =====\n');

for c = 1:nComp
    row      = runAllPairwiseTests(comparisons{c,1}, comparisons{c,2}, comparisons{c,3});
    allStats = [allStats; row];
end

% Bonferroni correction across 6 comparisons (applied per test family)
rawCols = {'MWU_P','Welch_P','Ftest_P','BF_P','Levene_P','KS_P','AD_P'};
adjCols = {'MWU_Adj','Welch_Adj','Ftest_Adj','BF_Adj','Levene_Adj','KS_Adj','AD_Adj'};
for t = 1:numel(rawCols)
    allStats.(adjCols{t}) = min(allStats.(rawCols{t}) * nComp, 1);
end

% Expose MWU p-values for significance brackets in Step 9
p_Bone_CR    = allStats.MWU_P(strcmp(allStats.Comparison, 'Bone Clone: TC vs Random'));
p_Bone_CA    = allStats.MWU_P(strcmp(allStats.Comparison, 'Bone Clone: TC vs Aligned'));
p_Bone_RA    = allStats.MWU_P(strcmp(allStats.Comparison, 'Bone Clone: Random vs Aligned'));
p_Parental_CR = allStats.MWU_P(strcmp(allStats.Comparison, 'Parental: TC vs Random'));
p_Parental_CA = allStats.MWU_P(strcmp(allStats.Comparison, 'Parental: TC vs Aligned'));
p_Parental_RA = allStats.MWU_P(strcmp(allStats.Comparison, 'Parental: Random vs Aligned'));

% ------------------ Step 8: Save All Statistical Results ------------------

pairwiseFile  = fullfile(mainFolder, 'BioRep1_Comprehensive_Stats_CellLevel.csv');
normalityFile = fullfile(mainFolder, 'BioRep1_Normality_Results.csv');

writetable(allStats,        pairwiseFile);
writetable(normalityResults, normalityFile);

fprintf('\nPairwise statistics saved to:\n  %s\n', pairwiseFile);
fprintf('Normality results saved to:\n  %s\n',     normalityFile);
fprintf('\nAnalysis complete.\n');

% =========================================================================
% Step 9: Add Significance Brackets to Box Plot and Violin Plot
% =========================================================================
% Each row of bracketData = { x_left,  x_right,  p_value }
% x positions follow the group order:
%   1 = Bone Clone-TC      4 = Parental-TC
%   2 = Bone Clone-Random  5 = Parental-Random
%   3 = Bone Clone-Aligned 6 = Parental-Aligned

bracketData = { ...
    1, 2, p_Bone_CR;      % Bone Clone:  TC vs Random
    2, 3, p_Bone_RA;      % Bone Clone:  Random vs Aligned
    1, 3, p_Bone_CA;      % Bone Clone:  TC vs Aligned
    4, 5, p_Parental_CR;  % Parental:    TC vs Random
    5, 6, p_Parental_RA;  % Parental:    Random vs Aligned
    4, 6, p_Parental_CA;  % Parental:    TC vs Aligned
};

figNames = {'Boxplot - Aspect Ratio', 'Violin Plot - Aspect Ratio'};
for f = 1:numel(figNames)
    hFig = findobj('Type', 'figure', 'Name', figNames{f});
    if isempty(hFig)
        fprintf('Figure "%s" not found — skipping brackets.\n', figNames{f});
        continue;
    end
    figure(hFig);
    ax = gca;
    hold(ax, 'on');
    drawSignificanceBrackets(ax, bracketData);
    hold(ax, 'off');
end

% =========================================================================
% Local Functions
% =========================================================================

% ------------------ Original function (unchanged) ------------------
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
% =========================================================================
% New Local Functions — Comprehensive Statistical Tests
% =========================================================================

% ------------------------------------------------------------------
function results = runAllPairwiseTests(g1, g2, label)
% RUNALLPAIRWISETESTS  Seven tests across four families on one pair of groups.
%
%   CENTRE  → Mann-Whitney U | Welch t-test (log-transformed)
%   SPREAD  → F-test         | Brown-Forsythe  | Levene's
%   SHAPE   → KS test        | Anderson-Darling (2-sample)

    g1 = g1(~isnan(g1) & g1 > 0);
    g2 = g2(~isnan(g2) & g2 > 0);

    fprintf('\n--- %s  (n1 = %d | n2 = %d) ---\n', label, numel(g1), numel(g2));

    % ---- CENTRE: location / mean ----------------------------------------

    % 1. Mann-Whitney U (ranksum) — non-parametric median comparison
    if numel(g1) < 2 || numel(g2) < 2
        p_mwu = NaN;
    else
        [p_mwu, ~] = ranksum(g1, g2);
    end
    fprintf('  [Centre]  Mann-Whitney U          p = %.5f\n', p_mwu);

    % 2. Welch t-test on log-transformed data — parametric mean comparison
    %    Uses log() because aspect ratios are right-skewed and multiplicative;
    %    log-transform often achieves approximate normality.
    log1 = log(g1);   log2 = log(g2);
    if numel(log1) < 2 || numel(log2) < 2
        p_welch = NaN;
    else
        [~, p_welch] = ttest2(log1, log2, 'Vartype', 'unequal');
    end
    fprintf('  [Centre]  Welch t-test (log)       p = %.5f\n', p_welch);

    % ---- SPREAD: variance -----------------------------------------------

    % 3. F-test (vartest2) — assumes normality; sensitive to outliers
    if numel(g1) < 2 || numel(g2) < 2
        p_f = NaN;
    else
        [~, p_f] = vartest2(g1, g2);
    end
    fprintf('  [Spread]  F-test (vartest2)        p = %.5f\n', p_f);

    % 4. Brown-Forsythe — median-based; robust to skew and outliers
    p_bf = brownForsytheTest(g1, g2);
    fprintf('  [Spread]  Brown-Forsythe           p = %.5f\n', p_bf);

    % 5. Levene's test — mean-based; standard homogeneity-of-variance test
    p_lev = levenesTest(g1, g2);
    fprintf('  [Spread]  Levene''s test            p = %.5f\n', p_lev);

    % ---- SHAPE: full distribution ----------------------------------------

    % 6. Two-sample KS test — sensitive to shifts in median/location
    if numel(g1) < 2 || numel(g2) < 2
        p_ks = NaN;
    else
        [~, p_ks] = kstest2(g1, g2);
    end
    fprintf('  [Shape]   KS test                  p = %.5f\n', p_ks);

    % 7. Two-sample Anderson-Darling — more sensitive to tail differences
    %    than KS; useful for detecting rare super-elongated cells.
    [p_ad, ~] = adtest2Sample(g1, g2);
    fprintf('  [Shape]   Anderson-Darling         p = %.5f\n', p_ad);

    results = table( ...
        {label}, p_mwu, p_welch, p_f, p_bf, p_lev, p_ks, p_ad, ...
        'VariableNames', { ...
            'Comparison', ...
            'MWU_P',   'Welch_P', ...
            'Ftest_P', 'BF_P', 'Levene_P', ...
            'KS_P',    'AD_P' });
end

% ------------------------------------------------------------------
function p = brownForsytheTest(g1, g2)
% BROWNFORSYTHETEST  Median-based Levene test (Brown-Forsythe, 1974).
%   Transforms each observation to its absolute deviation from the group
%   median, then compares the two deviation sets with a standard t-test.
%   Robust to non-normality and outliers — preferred for skewed data.

    g1 = g1(:);   g2 = g2(:);
    z1 = abs(g1 - median(g1));
    z2 = abs(g2 - median(g2));
    if numel(z1) < 2 || numel(z2) < 2
        p = NaN;  return;
    end
    [~, p] = ttest2(z1, z2);   % equal-variance 2nd-stage t-test (standard)
end

% ------------------------------------------------------------------
function p = levenesTest(g1, g2)
% LEVENESTEST  Mean-based Levene test (Levene, 1960).
%   Transforms each observation to its absolute deviation from the group
%   mean, then compares the two deviation sets with a standard t-test.

    g1 = g1(:);   g2 = g2(:);
    z1 = abs(g1 - mean(g1));
    z2 = abs(g2 - mean(g2));
    if numel(z1) < 2 || numel(z2) < 2
        p = NaN;  return;
    end
    [~, p] = ttest2(z1, z2);   % equal-variance 2nd-stage t-test (standard)
end

% ------------------------------------------------------------------
function [p, ADstat] = adtest2Sample(x, y)
% ADTEST2SAMPLE  Two-sample Anderson-Darling test.
%   Based on Scholz & Stephens (1987).  More sensitive than the KS test
%   in the distribution tails — useful for detecting rare super-elongated
%   or near-circular subpopulations.
%
%   P-value via chi-squared(1) asymptotic approximation; reliable for
%   n > ~20 per group, which is always the case for cell-level data.

    x = x(:);  y = y(:);
    m = numel(x);  n = numel(y);  N = m + n;

    if m < 2 || n < 2
        p = NaN;  ADstat = NaN;  return;
    end

    Z = sort([x; y]);   % combined order statistics

    ADstat = 0;
    for j = 1:N-1
        fj = j / N;                     % pooled ECDF value at Z(j)
        if fj <= 0 || fj >= 1, continue; end
        Fx = sum(x <= Z(j)) / m;        % ECDF of x at Z(j)
        Fy = sum(y <= Z(j)) / n;        % ECDF of y at Z(j)
        ADstat = ADstat + (Fx - Fy)^2 / (fj * (1 - fj));
    end
    ADstat = (m * n / N) * ADstat;

    % Asymptotic p-value: chi-squared with 1 df (k = 2 samples → df = k-1)
    p = 1 - chi2cdf(ADstat, 1);
    p = max(1e-6, min(p, 1));   % clamp to valid probability range
end
% ------------------ New: draw all brackets on an axes ------------------
function drawSignificanceBrackets(ax, bracketData)
% DRAWSIGNIFICANCEBRACKETS  Overlays significance brackets on a MATLAB axes.
%
%   bracketData  – n×3 cell array: {x1, x2, p_value} per row.
%
%   Brackets are stacked using an interval-graph colouring algorithm so
%   horizontally overlapping spans never share the same height level.
%   The y-axis is expanded automatically if brackets exceed the current limit.

    % ---- Sort by span (shorter brackets at lower levels for cleaner layout)
    n     = size(bracketData, 1);
    spans = zeros(n, 1);
    for i = 1:n
        spans(i) = bracketData{i,2} - bracketData{i,1};
    end
    [~, sortIdx] = sort(spans, 'ascend');
    bracketData  = bracketData(sortIdx, :);

    % ---- Assign a stacking level to each bracket (greedy interval colouring)
    levels = assignBracketLevels(bracketData);

    % ---- Geometry (relative to current y-axis)
    % ---- Detect background colour and pick a contrasting bracket colour
    bgCol  = ax.Color;                          % axes background
    isDark = mean(bgCol) < 0.5;                 % true for dark themes
    if isDark
        bracketCol = [1.00 1.00 1.00];          % white on dark background
    else
        bracketCol = [0.15 0.15 0.15];          % near-black on light background
    end

    yLim      = ax.YLim;
    yRange    = yLim(2) - yLim(1);
    yBase     = yLim(2) - yRange * 0.18;   % lowest bracket sits here
    stepSize  = yRange * 0.08;             % vertical gap between levels
    tickH     = yRange * 0.025;            % height of the L-shaped tick
    textGap   = yRange * 0.010;            % gap between bracket top and star text

    % Expand y-axis if needed to fit all brackets + text
    maxLevel  = max(levels);
    yRequired = yBase + maxLevel * stepSize + tickH + textGap + yRange * 0.06;
    if yRequired > yLim(2)
        ylim(ax, [yLim(1), yRequired]);
    end

    % ---- Draw each bracket
    for i = 1:n
        x1 = bracketData{i,1};
        x2 = bracketData{i,2};
        p  = bracketData{i,3};

        if isnan(p), continue; end

        starStr = pToStars(p);
        lvl     = levels(i);
        y0      = yBase + lvl * stepSize;        % foot of the bracket
        y1      = y0   + tickH;                  % top of the bracket

        % Bracket shape: two vertical ticks connected by a horizontal bar
        line(ax, [x1, x1, x2, x2], [y0, y1, y1, y0], ...
            'Color',     bracketCol, ...
            'LineWidth', 1.2, ...
            'Clipping',  'off');

        % Significance label centred above the bracket
        text(ax, (x1 + x2) / 2, y1 + textGap, starStr, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize',   9, ...
            'FontWeight', 'bold', ...
            'Color',      bracketCol);
    end
end

% ------------------ New: greedy level assignment ------------------
function levels = assignBracketLevels(bracketData)
% Assign each bracket the lowest non-negative integer level such that no
% two horizontally overlapping brackets share the same level.
% A small tolerance (0.1 x-units) prevents brackets sharing an endpoint
% from being counted as overlapping.

    n      = size(bracketData, 1);
    levels = zeros(n, 1);

    for i = 1:n
        x1i = bracketData{i,1};
        x2i = bracketData{i,2};

        % Collect levels already occupied by overlapping brackets
        occupied = [];
        for j = 1:i-1
            x1j = bracketData{j,1};
            x2j = bracketData{j,2};
            overlaps = (x1i < x2j - 0.1) && (x2i > x1j + 0.1);
            if overlaps
                occupied(end+1) = levels(j); %#ok<AGROW>
            end
        end

        % Assign lowest free level
        lvl = 0;
        while ismember(lvl, occupied)
            lvl = lvl + 1;
        end
        levels(i) = lvl;
    end
end

% ------------------ New: p-value → significance stars ------------------
function s = pToStars(p)
    if     p < 0.0001,  s = '****';
    elseif p < 0.001,   s = '***';
    elseif p < 0.01,    s = '**';
    elseif p < 0.05,    s = '*';
    else,               s = 'ns';
    end
end
% Save your data to the main project folder
excelFileName = fullfile(mainFolder, 'BioRep1_CellData.xlsx');
writetable(allCellData, excelFileName);
fprintf('Success! Data saved to Excel for Python analysis.\n');