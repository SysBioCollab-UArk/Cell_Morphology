%% TIBD Morphological Analysis - Bone Clone (Unaligned)
clear, clc, close all

% ------------------ Step 1: Configuration ------------------
main_folder = '/Users/hilmerdiaz/PycharmProjects/Cell_Morphology/AspectRatio/BioRep2_24Hrs_500kP_300kB_Feb_2026';
outDir = fullfile(main_folder, 'Cluster_Results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cp = cellpose(model="cyto2");
allAngles = []; 

% ------------------ Step 2: Nested Cluster Loop (Batch of 5x5) ------------------
for group = 8:9
    for imgNum = 1:5
        fileName = sprintf('%d_%d.tif', group, imgNum);
        fullPath = fullfile(main_folder, fileName);
        
        if ~exist(fullPath, 'file')
            fprintf('Skipping %s (file not found)\n', fileName);
            continue;
        end
        
        fprintf('Processing Cluster Image: %s\n', fileName);
        
        % Load and Pre-process
        img_raw = imread(fullPath);
        img = mean(single(img_raw), 3);
        img_adj = adapthisteq(img./65535);
        
        % Segment
        cellLabel = segmentCells2D(cp, img_adj, ImageCellDiameter=40);
        props = regionprops(cellLabel, img, 'Area', 'MajorAxisLength', ...
            'MinorAxisLength', 'Orientation', 'Circularity', 'Centroid');
        
        if isempty(props), continue; end
        
        % ------------------ Step 3: Filtering & Numbering ------------------
        % Create a base for the numbered TIF image
        labelImg = uint8((cellLabel > 0)) .* 255;
        labelImg = repmat(labelImg, 1, 1, 3); 

        for i = 1:length(props)
            props(i).Number = i; % Synchronize with CSV
            props(i).AspectRatio = props(i).MajorAxisLength / props(i).MinorAxisLength;
            
            % Apply the Isolation Thresholds
            if props(i).AspectRatio > 1.5 && props(i).Circularity < 0.6
                props(i).Class = 'Elongated_Bone_Clone';
                allAngles = [allAngles; deg2rad(props(i).Orientation)];
            elseif props(i).AspectRatio > 1.1
                props(i).Class = 'Alive_Rounded';
            else
                props(i).Class = 'Dead/Debris';
            end

            % DRAW THE NUMBER on the image
            labelImg = insertText(labelImg, round(props(i).Centroid), int2str(i), ...
                'AnchorPoint', 'Center', 'FontColor', 'w', 'BoxOpacity', 0, 'FontSize', 22);
        end
        
       % ------------------ Step 4: Color-Coded Numbering ------------------
        % Use the adjusted grayscale image as the background for the labels
        labelImg = uint8(img_adj .* 255); 
        labelImg = repmat(labelImg, 1, 1, 3); % Convert to RGB to allow for Green/Red text
        
        for i = 1:length(props)
            pos = props(i).Centroid;
            cellID = num2str(i);
            
            % Determine Color based on the Classification we did in Step 3
            if strcmp(props(i).Class, 'Elongated_Bone_Clone')
                textColor = 'green'; % These are the ones in your Rose Plot
            else
                textColor = 'red';   % These were excluded
            end
            
            % DRAW THE COLORED NUMBER
            labelImg = insertText(labelImg, pos, cellID, ...
                'AnchorPoint', 'Center', ...
                'FontSize', 28, ... % Slightly larger for visibility
                'TextColor', textColor, ...
                'BoxOpacity', 0);
        end
        
       imwrite(labelImg, fullfile(outDir, [fileName(1:end-4), '_labeledCells.tif']));
        
    end 
end 

% ------------------ Step 5: Overall Orientation & Index ------------------
if ~isempty(allAngles)
    S = abs(mean(exp(2i * allAngles))); 
    
    figure;
    polarhistogram(allAngles, 36, 'FaceColor', 'b', 'FaceAlpha', 0.5);
    title(['Orientation: Bone Clone Unaligned (Index: ', num2str(S, '%.2f'), ')']);
    subtitle(['Total Cells (N) across 25 images: ', num2str(length(allAngles))]);
    
    saveas(gcf, fullfile(outDir, 'Overall_Orientation_RosePlot.png'));
    
    fprintf('\n--- Final Cluster Stats ---\n');
    fprintf('Total Cells Picked Up: %d\n', length(allAngles));
    fprintf('Alignment Index (S): %.4f\n', S);
else
    disp('No data available for Index calculation.');
end