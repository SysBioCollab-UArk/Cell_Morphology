%% Cell Classification and Labeling Script
% This script processes TIFF images, segments cells, measures properties,
% classifies them as Alive or Dead, and overlays colored labels
% accordingly.

clear, clc, close all

% ------------------ Step 1: Set main folder ------------------
%mainfolder = '/Users/alexandravega/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/JUNE_JULY_2025/2025_07_10_09_52_25--24_Hrs_07_10_2025_ALL/1_second';
main_folder = '/Users/alexandravega/PycharmProjects/Cell_Morphology/AspectRatio';
%main_folder = '/Users/alexandravega/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/JUNE_JULY_2025/2025_07_10_09_52_25--24_Hrs_07_10_2025_ALL/June_July_Biorep2';
%main_folder = '/Users/alexandravega/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/MAY 2025/05_23_2025/MayImages_20x';
% ------------------ Step 2: Get all TIF images ------------------
tif_files = dir(fullfile(main_folder, '*.tif'));
%tif_files = dir(fullfile(main_folder, '*.tif'));
%names = {tif_files.name};

% first integer 1–5, then '_' int '_' int, ending .tif
%keep = ~cellfun('isempty', regexp(names, '^[1-5]_\d+_\d+\.tif$', 'once'));

% tif_files = tif_files(keep); % Filter tif_files based on the keep condition
% %for k = 1:numel(tif_files)
% %    fprintf('%s\n', tif_files(k).name);
% %end
% if isempty(tif_files)
%     disp('No matching .tif files found.');
% end

% ------------------ Step 3: Load cellpose model ------------------
cp = cellpose(model="cyto2");

% ------------------ Step 4: Loop through images ------------------
%densityTable = struct('Name',{},'Density',{},'umPerPixel',{});
%densityTable = struct('Name',{}, 'Density',{},'umPixel'{},'AvgIntensity'{};);
%densityTable = struct('Name',{},'Density',{},'umPerPixel'{},'AvgIntensity'{});
% Store density and average intensity in the density table
    %densityTable(k).Name = tif_files(k).name;
    %densityTable(k).Density = cellDensity;
    %densityTable(k).umPerPixel = microns_per_pixel;
    %densityTable(k).AvgIntensity = mean([props.MeanIntensity]);

%densityTable = struct('Name',{},'Density',{},'umPerPixel',{},'AvgIntensity',{});
%...
%densityTable(k).AvgIntensity = mean(img(:));


%microns_per_pixel = 0.5357154; % Calibration factor for Image

%for k = 1:length(tif_files)

 for k =1: min(5, length(tif_files))
     fprintf('Processing image %d of %d: %s\n', ...
        k, min(5,length(tif_files)), tif_files(k).name);


% Load an image
    img = imread(fullfile(tif_files(k).folder, tif_files(k).name));
    img = mean(single(img),3);
    %rawImg =mean(single(imread(fullfile(tif_files(k).folder, tif_files(k).name)))3); % added for intensity
    img = adapthisteq(img./65535);

    info = imfinfo(fullfile(tif_files(k).folder, tif_files(k).name));
    if ~strcmpi(info.ResolutionUnit, "Centimeter")
        a
    else
        microns_per_pixel = 1/(info.XResolution / 10000);
    end

    % Segment the image using cellpose
    cellLabel = segmentCells2D(cp, img, ImageCellDiameter=40);

    % Get cell properties (via regionprops)
    %props = regionprops(cellLabel, 'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength','Eccentricity');
    props = regionprops(cellLabel, img,'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', ...
        'MeanIntensity', 'MaxIntensity', 'MinIntensity'); % added for intensity


    %props = orderfields(props, {'Class','Number','CentroidX','CentroidY','Area','ActualArea','AspectRatio', ...
    %'MeanIntensity','MaxIntensity','MinIntensity', ...
    %'MajorAxisLength','MinorAxisLength','MajorAxis_um','MinorAxis_um','Eccentricity'});

    % Calculate aspect ratio and do some formatting
    aliveCount = 0;
    for i = 1:length(props)
        props(i).Number = i;
        props(i).AspectRatio = props(i).MajorAxisLength ./ props(i).MinorAxisLength;
        props(i).MajorAxis_um = props(i).MajorAxisLength * microns_per_pixel;
        props(i).MinorAxis_um = props(i).MinorAxisLength * microns_per_pixel;
        %props = regionprops(cellLabel, rawImg, 'MeanIntensity', 'MaxIntensity', 'MinIntensity'); % added for intensity

        if props(i).AspectRatio > 1.2
            props(i).Class = 'Alive'; 
            aliveCount = aliveCount + 1;

        else 
            props(i).Class = 'Dead'; 
        end

        props(i).CentroidX = props(i).Centroid(1);
        props(i).CentroidY = props(i).Centroid(2);
        props(i).ActualArea = props(i).Area*(microns_per_pixel^2); % px^2 -> um^2 % I added the ^2 inside the parenthesis to the micros_per_pixel
      
    end

    %Label and Visualize Cells
    labels = cellLabel;
    labelImg = uint8((labels > 0)).* 255; % Create binary mask for labeled regions
    labelImg = repmat(labelImg, 1,1, 3); % Convert to RGB


    for i = 1:length(props)
        labelImg = insertText(labelImg, round(props(i).Centroid), int2str(i),'AnchorPoint', 'Center', 'FontColor', 'r', 'BoxOpacity', 0, 'FontSize', 32); 
    end

    figure; imagesc(labelImg);
    title(['Labeled Cells - ', tif_files(k).name]);

    % Save output images
    outDir = fullfile(tif_files(k).folder, 'output');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    imwrite(labelImg, fullfile(outDir, [tif_files(k).name(1:end-4), '_labeledCells.tif']));

    %imwrite(uint16(cellLabel), fullfile(outDir, [tif_files(k).name(1:end-4), '_cellLabel.tif']));
    %imwrite(uint16(img .* 65535), fullfile(outDir, [tif_files(k).name(1:end-4), '_intAdj.tif']));

    % Write cell 
    % property data
    %writecell([fieldnames(props)'; struct2cell(props)], fullfile(outDir, [tif_files(k).name(1:end-4), '_cellProps.csv']));
    
    % Re-order fields for readability
    props = orderfields(props, {'Number','Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','MeanIntensity','MinIntensity','MaxIntensity','AspectRatio','MajorAxis_um','MinorAxis_um','Class','CentroidX','CentroidY','ActualArea'});
    a
    % --- Make sure props exists before ordering ---
    if ~exist('props','var') || isempty(props)
    warning('No cells found (props empty) for image: %s', tif_files(k).name);
    continue;
    end

    % Organize property structure
    props = rmfield(props,"Centroid");
    % Organize property structure
    props = rmfield(props, "Centroid");% added
    %props = orderfields(props, {...});
    props = orderfields(props, {'Class','Number','CentroidX','CentroidY','Area','ActualArea','AspectRatio','MajorAxisLength','MinorAxisLength','MajorAxis_um', 'MinorAxis_um', 'Eccentricity', 'MeanIntensity','MaxIntensity', 'MinIntensity'});
    props = regionprops(cellLabel, img, ...
    'Centroid','Area','MajorAxisLength','MinorAxisLength','Eccentricity', ...
    'MeanIntensity','MaxIntensity','MinIntensity'); % Added for intensity
    % Generate map of alive/dead cells
    aliveMap = zeros(size(img),'uint8');

    for i = 1:length(props)
        if strcmp(props(i).Class,'Alive')
            aliveMap(cellLabel==i) = 1;
        else
            aliveMap(cellLabel==i) = 2;
        end
    end

    % Calculate cell density
    cellDensity = (sum(logical(cellLabel(:))) ./ numel(cellLabel)) .* 100;

    % Save output images
    imwrite(uint16(cellLabel), fullfile(outDir, [tif_files(k).name(1:end-4), '_cellLabel.tif']));
    imwrite(uint16(img.*65535), fullfile(outDir,[tif_files(k).name(1:end-4),'_intAdj.tif']));
    %outDir = fullfile(tif_files(k).folder), 'output');
    %if ~exist(outDir,'dir')
        %mkdir(outDir)
    % Write Cell property data
    writecell([fieldnames(props)'; squeeze(struct2cell(props))'], fullfile(outDir, [tif_files(k).name(1:end-4),'_cellProps.csv']));

    % Add density to table
    densityTable(k).Name = tif_files(k).name;
    densityTable(k).Density =cellDensity;
    densityTable(k).umPerPixel = microns_per_pixel;
    densityTable(k).AvgIntensity = mean(img(:)); % Added for intensity


end

    % Write Summary density file
    %writecell([fieldnames(densityTable)'; struct2cell(densityTable)]', fullfile(tif_files(1).folder, 'output', 'cellDensity.scv'));
    writecell([fieldnames(densityTable)'; squeeze(struct2cell(densityTable))'], fullfile(tif_files(1).folder, 'output', 'cellDensity.csv'));
    
    
    %writecell([fieldnames(densityTable)'; squeeze(struct2cell(densityTable))], fullfile(tif_files(1).folder, 'output', 'cellDensity.csv'));
    

    %imwrite(uint16(cellLabel), fullfile, (outDir, [tif_files(k).folder,'/output/', tif_files(k).name(1:end-4),'_cellLabel.tif']);
    %imwrite(uint16(img.*65535), fullfile,(outDir), [tif_files(k).folder,'/output/', tif_files(k).name(1:end-4),'_intAdj.tif']);

    % Write cell property data
    %writecell([fieldnames(props)';squeeze(struct2cell(props))'], fullfile(outDir, [tif_files(k).folder,'/output/', tif_files(k).name(1:end-4),'_cellProps.csv']);

    % Add density to table
    %densityTable(k).Name = tif_files(k).name;
    %densityTable(k).Density = cellDensity;
%end

%writecell([fieldnames(densityTable)';squeeze(struct2cell(densityTable))'], fullfile([tif_files(k).folder,'output', 'cellDensity.csv']);