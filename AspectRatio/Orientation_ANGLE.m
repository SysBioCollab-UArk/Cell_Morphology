%% =========================================================================
%  TIBD Orientation Analysis Pipeline (Publication-Optimized) k
%  
%  Proves scaffold-guided cell alignment using orientation angle data.
% =========================================================================
clear; clc; close all;

%% ─────────────────────────────────────────────────────────────────────────
%  USER CONFIGURATION
% ─────────────────────────────────────────────────────────────────────────
%CSV_PATH = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/2026_04_09_10_38_48--BioRep2_Double_Cells_04_09_2026_Final/BioRep2_HighPrecision_ALL_DATA.csv';
CSV_PATH = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026/BioRep1_HighPrecision_ALL_DATA.csv';
%CSV_PATH = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep1_Double_03_24_2026/BioRep1_RAW_ALL_DATA.csv';   
OUT_DIR         = fullfile(fileparts(CSV_PATH), 'Orientation_Results');
ALPHA           = 0.05;
CONDITION_ORDER = {'TC', 'Random', 'Aligned'};
N_BINS          = 36;   % Rose plot bins (36 = every 5 degrees)

% Colors: Bone Clone = dark blue | Parental = orange
CL_COLORS = [0.15 0.35 0.70;
             0.85 0.45 0.10];

%% ─────────────────────────────────────────────────────────────────────────
%  LOAD & CLEAN
% ─────────────────────────────────────────────────────────────────────────
fprintf('=== Loading CSV ===\n');
T = readtable(CSV_PATH, 'TextType', 'string');

% Fix known capitalisation issue
T.Condition = replace(T.Condition, 'Tc', 'TC');

% Clean hidden characters
cleanStr    = @(col) regexprep(strtrim(col), '[^\x20-\x7E]', '');
T.CellLine  = cleanStr(T.CellLine);
T.Condition = cleanStr(T.Condition);

if ~ismember('Orientation', T.Properties.VariableNames)
    error('Orientation column not found in CSV.');
end

fprintf('Rows loaded: %d\n', height(T));
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

cellLines  = ["Bone Clone"; "Parental"];
conditions = string(CONDITION_ORDER);
nCL        = numel(cellLines);
nCond      = numel(conditions);

%% ─────────────────────────────────────────────────────────────────────────
%  STEP 1 — Convert Orientation for Circular Statistics
% ─────────────────────────────────────────────────────────────────────────
T.OrientDeg = T.Orientation;
T.OrientRad_doubled = deg2rad(2 .* T.OrientDeg);

%% ─────────────────────────────────────────────────────────────────────────
%  STEP 2 — Rose Plots
% ─────────────────────────────────────────────────────────────────────────
fprintf('=== STEP 2: Rose Plots ===\n');
fig_rose = figure('Name','Orientation Rose Plots','Color','w', 'Position',[50 50 1400 700]);

for cl = 1:nCL
    for co = 1:nCond
        spIdx = (cl-1)*nCond + co;
        ax    = subplot(nCL, nCond, spIdx);
        mask = T.CellLine == cellLines(cl) & T.Condition == conditions(co);
        vals = T.OrientDeg(mask);
        
        if numel(vals) < 5
            title(sprintf('%s\n%s\n(no data)', cellLines(cl), conditions(co)));
            continue;
        end
        
        valsRad    = deg2rad(vals + 0);          
        valsRadFull = [valsRad; valsRad + pi];    
        
        pax = polaraxes('Position', ax.Position);
        delete(ax);
        
        h = polarhistogram(valsRadFull, N_BINS*2, 'Normalization','probability', ...
            'FaceColor', CL_COLORS(cl,:), 'EdgeColor', 'w', 'FaceAlpha', 0.85, 'LineWidth', 0.5);
            
        r_val  = circ_r_axial(vals);
        mu_deg = circ_mean_axial(vals);
        mu_rad = deg2rad(mu_deg + 90);
        
        hold(pax,'on');
        maxR = max(h.Values);
        % Draw bold mean direction line extending slightly past the data
        polarplot(pax, [mu_rad mu_rad+pi], [maxR*1.15 maxR*1.15], '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 2.5);
        hold(pax,'off');
        
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        pax.ThetaTick = [0 45 90 135 180 225 270 315];
        pax.ThetaTickLabel = {'-90°','-45°','0°','45°','90°','135°','180°','-135°'};
        
        title(pax, sprintf('%s | %s\nn=%d  r=%.3f  μ=%.1f°', cellLines(cl), conditions(co), sum(mask), r_val, mu_deg), 'FontSize', 10, 'FontWeight','bold');
    end
end
sgtitle('Cell Orientation Rose Plots', 'FontSize', 16, 'FontWeight','bold');
exportgraphics(fig_rose, fullfile(OUT_DIR,'Rose_Plots.png'), 'Resolution',300);

%% ─────────────────────────────────────────────────────────────────────────
%  STEP 3 — Circular Statistics per Group
% ─────────────────────────────────────────────────────────────────────────
fprintf('=== STEP 3: Circular Statistics ===\n');
fprintf('%-14s %-10s %6s %8s %8s %8s %10s %12s\n', ...
    'CellLine','Condition','N','Mean_Ang','r','AI (%)','CircSD','Rayleigh_p');
fprintf('%s\n', repmat('-',1,85));
statsRows = {};

for cl = 1:nCL
    for co = 1:nCond
        mask = T.CellLine == cellLines(cl) & T.Condition == conditions(co);
        vals = T.OrientDeg(mask);
        n    = numel(vals);
        if n < 5, continue; end
        
        mu    = circ_mean_axial(vals);     
        r     = circ_r_axial(vals);        
        csd   = circ_std_axial(vals);      
        rp    = rayleigh_test_axial(vals); 
        aiPct = r * 100;                   
        
        % Prints the Alignment Index (aiPct) to the command window
        fprintf('%-14s %-10s %6d %7.2f° %8.4f %7.1f%% %7.2f° %12.4e\n', ...
            cellLines(cl), conditions(co), n, mu, r, aiPct, csd, rp);
            
        statsRows(end+1,:) = {char(cellLines(cl)), char(conditions(co)), ...
            n, mu, r, csd, rp, aiPct, ternary_str(rp < ALPHA, 'Significant','Not significant')}; %#ok<SAGROW>
    end
end

circStatsTable = cell2table(statsRows, 'VariableNames', ...
    {'CellLine','Condition','N','MeanAngle_deg','ResultantLength_r', ...
     'CircularSD_deg','Rayleigh_p','AlignmentIndex_pct','Alignment'});
writetable(circStatsTable, fullfile(OUT_DIR,'Circular_Statistics.csv'));

%% ─────────────────────────────────────────────────────────────────────────
%  STEP 4 — Watson U² Test (Robust Implementation)
% ─────────────────────────────────────────────────────────────────────────
watsonRows = {};
condPairs  = nchoosek(1:nCond, 2);

for cl = 1:nCL
    for p = 1:size(condPairs,1)
        co1 = condPairs(p,1); co2 = condPairs(p,2);
        mask1 = T.CellLine==cellLines(cl) & T.Condition==conditions(co1);
        mask2 = T.CellLine==cellLines(cl) & T.Condition==conditions(co2);
        vals1 = T.OrientDeg(mask1); vals2 = T.OrientDeg(mask2);
        
        if numel(vals1)<5 || numel(vals2)<5, continue; end
        
        U2  = watson_u2_axial(vals1, vals2);
        sig = watson_interpret(U2);
        watsonRows(end+1,:) = {char(cellLines(cl)), char(conditions(co1)), char(conditions(co2)), U2, sig}; %#ok<SAGROW>
    end
end
watsonTable = cell2table(watsonRows, 'VariableNames', {'CellLine','Condition_1','Condition_2','Watson_U2','Interpretation'});
writetable(watsonTable, fullfile(OUT_DIR,'Watson_U2_Tests.csv'));

%% ─────────────────────────────────────────────────────────────────────────
%  STEP 5 — Alignment Index Bar Chart (Dynamically Scaled)
% ─────────────────────────────────────────────────────────────────────────
fprintf('=== STEP 5: Alignment Index Plot ===\n');
fig_ai = figure('Name','Alignment Index','Color','w','Position',[100 100 800 550]);
hold on;
barWidth = 0.32; xBase = 1:nCond;
maxAI = 0; % Track maximum value for dynamic scaling

for cl = 1:nCL
    rVals = zeros(1,nCond);
    for co = 1:nCond
        row = strcmp(circStatsTable.CellLine, char(cellLines(cl))) & ...
              strcmp(circStatsTable.Condition, char(conditions(co)));
              
        % If matches exist, grab ONLY the first one to prevent array size errors
        if any(row)
            allMatches = circStatsTable.AlignmentIndex_pct(row);
            rVals(co) = allMatches(1); 
        end
    end
    
    maxAI = max(maxAI, max(rVals));
    
    xOff = (cl-1.5)*barWidth;
    b = bar(xBase+xOff, rVals, barWidth, 'FaceColor',CL_COLORS(cl,:),'EdgeColor','none','FaceAlpha',0.9,'DisplayName',char(cellLines(cl)));
    
    % Add subtle values above the bars
    for i = 1:length(rVals)
        text(xBase(i)+xOff, rVals(i) + (maxAI*0.02), sprintf('%.1f%%', rVals(i)), 'HorizontalAlignment','center','FontSize',9,'Color',[0.3 0.3 0.3]);
    end
end

set(gca,'XTick',xBase,'XTickLabel',conditions,'FontSize',12,'LineWidth',1);
ylabel('Alignment Index (%)', 'FontSize',13, 'FontWeight','bold');
title('Scaffold-Guided Alignment Index', 'FontSize',15,'FontWeight','bold');
legend('Location','northeast','FontSize',11);

% Dynamic Y-Axis Scaling (gives breathing room above the tallest bar)
ylim([0, maxAI * 1.25]);
grid on; box on; set(gca,'GridAlpha',0.2); hold off;
exportgraphics(fig_ai, fullfile(OUT_DIR,'Alignment_Index.png'),'Resolution',300);

%% ─────────────────────────────────────────────────────────────────────────
%  STEP 6 — Linear Density Distributions (with KDE Smoothing)
% ─────────────────────────────────────────────────────────────────────────
fprintf('=== STEP 6: Linear Smoothed Histograms ===\n');
fig_hist = figure('Name','Linear Orientation Histograms','Color','w', 'Position',[100 100 1400 450]);

% First pass: find the global maximum probability density to unify the Y-axes
globalMaxY = 0;
for co = 1:nCond
    for cl = 1:nCL
        mask = T.CellLine==cellLines(cl) & T.Condition==conditions(co);
        vals = T.OrientDeg(mask);
        if numel(vals) > 5
            [f, ~] = ksdensity(vals, linspace(-90,90,100));
            globalMaxY = max(globalMaxY, max(f));
        end
    end
end

for co = 1:nCond
    subplot(1,nCond,co); hold on;
    lgd = {};
    for cl = 1:nCL
        mask = T.CellLine==cellLines(cl) & T.Condition==conditions(co);
        vals = T.OrientDeg(mask);
        if numel(vals) < 5, continue; end
        
        % Faded background histogram
        histogram(vals, linspace(-90,90,N_BINS+1), 'Normalization','pdf', 'FaceColor',CL_COLORS(cl,:), 'EdgeColor','none', 'FaceAlpha',0.2);
        
        % Bold Kernel Density smoothing curve
        [f, xi] = ksdensity(vals, linspace(-90,90,200));
        plot(xi, f, 'Color', CL_COLORS(cl,:), 'LineWidth', 3.0);
        
        mu = circ_mean_axial(vals);
        lgd{end+1} = sprintf('%s (μ=%.1f°)', cellLines(cl), mu); %#ok<SAGROW>
    end
    xlabel('Orientation Angle (°)','FontSize',12);
    if co == 1, ylabel('Probability Density','FontSize',12); end
    title(sprintf('%s Scaffold',conditions(co)),'FontSize',14,'FontWeight','bold');
    
    % Use unified Y-axis for honest visual comparison
    ylim([0, globalMaxY * 1.1]);
    xlim([-90 90]); xticks(-90:30:90);
    
    % Format legends to only show the solid KDE lines
    hLines = findobj(gca, 'Type', 'Line');
    legend(flipud(hLines), lgd, 'Location','northwest','FontSize',10);
    grid on; box on; set(gca,'GridAlpha',0.2); hold off;
end
sgtitle('Cell Orientation Density Distributions','FontSize',16,'FontWeight','bold');
exportgraphics(fig_hist, fullfile(OUT_DIR,'Orientation_Histograms_Linear.png'),'Resolution',300);

%% =========================================================================
%  CIRCULAR STATISTICS HELPER FUNCTIONS (Robust & Corrected)
% =========================================================================
function mu = circ_mean_axial(deg)
    rad = deg2rad(2 .* deg);
    mu_rad = angle(mean(exp(1i .* rad))) / 2;
    mu = rad2deg(mu_rad);
end

function r = circ_r_axial(deg)
    rad = deg2rad(2 .* deg);
    r = abs(mean(exp(1i .* rad)));
end

function sd = circ_std_axial(deg)
    r = circ_r_axial(deg);
    sd = rad2deg(sqrt(-2 .* log(r))) / 2;
end

function p = rayleigh_test_axial(deg)
    n = numel(deg);
    r = circ_r_axial(deg);
    R = n .* r;
    Z = R.^2 ./ n;
    p = exp(-Z) .* (1 + (2.*Z - Z.^2)./(4.*n) - (24.*Z - 132.*Z.^2 + 76.*Z.^3 - 9.*Z.^4)./(288.*n.^2));
    p = max(0, min(1, p));
end

function U2 = watson_u2_axial(deg1, deg2)
    % SAFE TAGGING IMPLEMENTATION for massive N sizes
    a1 = mod(deg2rad(2.*deg1), 2*pi);
    a2 = mod(deg2rad(2.*deg2), 2*pi);
    n1 = numel(a1); n2 = numel(a2); n = n1+n2;
    
    taggedData = [a1(:), ones(n1,1); a2(:), zeros(n2,1)];
    sortedData = sortrows(taggedData, 1);
    
    isMember1 = sortedData(:, 2) == 1;
    C1 = cumsum(isMember1) ./ n1;
    C2 = cumsum(~isMember1) ./ n2;
    d  = C1 - C2;
    
    U2 = (n1.*n2./n^2) .* sum(d.^2 - mean(d).^2);
end

function s = watson_interpret(U2)
    if U2 > 0.267,      s = 'p < 0.01  (significant)';
    elseif U2 > 0.187,  s = 'p < 0.05  (significant)';
    else,               s = 'p > 0.05  (not significant)';
    end
end

function s = ternary_str(cond, a, b)
    if cond, s = a; else, s = b; end
end