clear, clc, close all

% ------------------ Step 1: Set main folder ------------------
%main_folder = '/Users/alexandravega/PycharmProjects/Cell_Morphology';
main_folder = '/Users/alexandravega/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/JUNE_JULY_2025/2025_07_10_09_52_25--24_Hrs_07_10_2025_ALL/June_July_Biorep2';
%main_folder = '/Users/alexandravega/Library/CloudStorage/Box-Box/UArk Sys Bio Collab/Projects/TIBD/qPCR/MAY 2025/05_23_2025/MayImages_20x';
% ------------------ Step 2: Get all TIF images ------------------
tif_files = dir(fullfile(main_folder, '*.tif'));
names = {tif_files.name};

% first integer 1–5, then '_' int '_' int, ending .tif
keep = ~cellfun('isempty', regexp(names, '^[1-5]_\d+_\d+\.tif$', 'once'));

tif_files = tif_files(keep); % Filter tif_files based on the keep condition
%for k = 1:numel(tif_files)
%    fprintf('%s\n', tif_files(k).name);
%end
if isempty(tif_files)
    disp('No matching .tif files found.');
end

% ------------------ Step 3: Load cellpose model ------------------
cp = cellpose(model="cyto2");
% ------------------ Step 4: Loop through images ------------------
densityTable = struct('Name',{},'Density',{}); 
for k = 1:length(tif_files)
 
% Load an image
    img = imread([tif_files(k).folder,'/',tif_files(k).name]);
    img = mean(single(img),3);
    img = adapthisteq(img./65535);

    % Segment the image using cellpose
    cellLabel = segmentCells2D(cp, img, ImageCellDiameter=40);

    % Get cell properties (via regionprops)
    props = regionprops(cellLabel, 'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength');

    % Calculate aspect ratio and do some formatting
    aliveCount = 0;
    for i = 1:length(props)
        props(i).Number = i;
        props(i).AspectRatio = props(i).MajorAxisLength ./ props(i).MinorAxisLength;
        if props(i).AspectRatio > 1.2
            props(i).Class = 'Alive'; 
            aliveCount = aliveCount + 1;

        else 
            props(i).Class = 'Dead'; 
        end

        props(i).CentroidX = props(i).Centroid(1);
        props(i).CentroidY = props(i).Centroid(2);
        props(i).ActualArea = props(i).Area*(0.5357154*0.5357154); % px^2 -> um^2
    end
    props = rmfield(props,"Centroid");
    props = orderfields(props, {'Class','Number','CentroidX','CentroidY','Area','ActualArea','AspectRatio','MajorAxisLength','MinorAxisLength'});

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
    if ~exist([tif_files(k).folder,'/output'],'dir')
        mkdir([tif_files(k).folder,'/output'])
    end

    imwrite(uint16(cellLabel), [tif_files(k).folder,'/output/', tif_files(k).name(1:end-4),'_cellLabel.tif']);
    imwrite(uint16(img.*65535), [tif_files(k).folder,'/output/', tif_files(k).name(1:end-4),'_intAdj.tif']);

    % Write cell property data
    writecell([fieldnames(props)';squeeze(struct2cell(props))'], [tif_files(k).folder,'/output/', tif_files(k).name(1:end-4),'_cellProps.csv']);

    % Add density to table
    densityTable(k).Name = tif_files(k).name;
    densityTable(k).Density = cellDensity;
end

writecell([fieldnames(densityTable)';squeeze(struct2cell(densityTable))'], [tif_files(k).folder,'/output/cellDensity.csv']);
