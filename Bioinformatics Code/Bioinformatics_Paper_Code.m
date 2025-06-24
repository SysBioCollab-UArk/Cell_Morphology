% Demo.m
%% define a parameter set
tiledim = 30;

lambda = 5; minSizeMSER = 30; maxSizeMSER = 4000; maxVariation = 1;

maxEcc = .7; minSizeSplit = 30; maxSizeSplit = 1000;


%load an image
I_org = imread('Demo1.png');
%I_org = uint8(I_org);

%run the code and visualize it
bw = segmentImage(I_org,'visualize',true);





% bgest.m
function [ZI,interp]=bgest(I,binsize)

%% init
fprintf('init...\n')

imgsize=size(I);

x = nan(floor((imgsize(1)-binsize)/binsize*2*(imgsize(2)-binsize)/binsize*2),1);
y = nan(floor((imgsize(1)-binsize)/binsize*2*(imgsize(2)-binsize)/binsize*2),1);
z = nan(floor((imgsize(1)-binsize)/binsize*2*(imgsize(2)-binsize)/binsize*2),1);
featuremat = nan(floor((imgsize(1)-binsize)/binsize*2*(imgsize(2)-binsize)/binsize*2),5);


imgwidth=size(I,2);
imghight=size(I,1);


counter=0;
fprintf('getting tiling features...\n')
for i = 1:binsize/2:imgsize(1)-binsize
    
    for j=1:binsize/2:imgsize(2)-binsize
        
        sub = I(i:i+binsize, j:j+binsize);
        counter=counter+1;
        y(counter) = i+binsize/2;
        x(counter) = j+binsize/2;
        z(counter)= mean(sub(:));
        featuremat(counter,:)=([std(sub(:)) skewness(sub(:)) max(sub(:))/min(sub(:)) kurtosis(sub(:)) var(sub(:))/mean(sub(:))  ]);
        
    end
end
%% cluster it
fprintf('clustering %d points...\n',counter)

[classes,type]=dbscan(featuremat,size(featuremat,2)+1,[]);

if numel(unique(classes))==1
    daclass = unique(classes);
else
    classstd = [];
    for c = unique(classes)
        if numel(featuremat(classes == c & type == 1,1))>200
            classstd(end+1)= mean(featuremat(classes == c,1));
        else
            classstd(end+1) = inf;
        end
    end
    [~,daclass] = min(classstd(2:end));
    daclass = daclass+1;
    daclasstemp  = unique(classes);
    daclass = daclasstemp(daclass);
end
interp=sum(classes == daclass& type == 1);

fprintf('using %d interpolation points...\n',interp)

%% interpolate

[XI YI] = meshgrid(1:imgwidth,1:imghight);
F=TriScatteredInterp(x(classes == daclass& type == 1),y(classes ==daclass& type == 1) ,z(classes ==daclass& type == 1),'natural');
ZI=F(XI,YI);

%% extrapolate
for i=find(sum(~isnan(ZI(1:imghight,:)))>1)
    ZI(:,i)=interp1(find(~isnan(ZI(:,i))),ZI(~isnan(ZI(:,i)),i), 1:imghight,'linear','extrap');
end
%%
for i=find(sum(~isnan(ZI(:,1:imgwidth))')>1)
    ZI(i,:)=interp1(find(~isnan(ZI(i,:))),ZI(i,~isnan(ZI(i,:))), 1:imgwidth,'linear','extrap');
end


%% fix strange extrapolations
ZI(ZI<min(I(:)))=min(I(:));
ZI(ZI>max(I(:)))=max(I(:));


fprintf('done\n')
end




% dbscan.m
function [class,type]=dbscan(x,k,Eps)

[m,n]=size(x);

if nargin<3 | isempty(Eps)
   [Eps]=epsilon(x,k);
end

x=[[1:m]' x];
[m,n]=size(x);
type=zeros(1,m);
no=1;
touched=zeros(m,1);

for i=1:m
    if touched(i)==0;
       ob=x(i,:);
       D=dist(ob(2:n),x(:,2:n));
       ind=find(D<=Eps);
    
       % more than 1 neighbour but less then k
       % thats noise
       % type outlier
       if length(ind)>1 & length(ind)<k+1       
          type(i)=0;
          class(i)=0;
       end
       % this is a single noise group?
       % type -1 outlier
       if length(ind)==1
          type(i)=-1;
          class(i)=-1;  
          touched(i)=1;
       end

       if length(ind)>=k+1; 
          type(i)=1;
          class(ind)=ones(length(ind),1)*max(no);
          
          while ~isempty(ind)
                ob=x(ind(1),:);
                touched(ind(1))=1;
                ind(1)=[];
                D=dist(ob(2:n),x(:,2:n));
                i1=find(D<=Eps);
     
                if length(i1)>1
                   class(i1)=no;
                   if length(i1)>=k+1;
                      type(ob(1))=1;
                   else
                      type(ob(1))=0;
                   end

                   for i=1:length(i1)
                       if touched(i1(i))==0
                          touched(i1(i))=1;
                          ind=[ind i1(i)];   
                          class(i1(i))=no;
                       end                    
                   end
                end
          end
          no=no+1; 
       end
   end
end

i1=find(class==0);
class(i1)=-1;
type(i1)=-1;
end

%...........................................
function [Eps]=epsilon(x,k)

% Function: [Eps]=epsilon(x,k)
%
% Aim: 
% Analytical way of estimating neighborhood radius for DBSCAN
%
% Input: 
% x - data matrix (m,n); m-objects, n-variables
% k - number of objects in a neighborhood of an object
% (minimal number of objects considered as a cluster)



[m,n]=size(x);

Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);
end

%............................................
function [D]=dist(i,x)

% function: [D]=dist(i,x)
%
% Aim: 
% Calculates the Euclidean distances between the i-th object and all objects in x	 
%								    
% Input: 
% i - an object (1,n)
% x - data matrix (m,n); m-objects, n-variables	    
%                                                                 
% Output: 
% D - Euclidean distance (m,1)



[m,n]=size(x);
D=sqrt(sum((((ones(m,1)*i)-x).^2)'));

if n==1
   D=abs((ones(m,1)*i-x))';
end
end




% segmentImage.m
%% load an image that was acquired in out-of-focus settings
function bw2 = segmentImage(I_org,varargin)


%%parse the parameters
p = inputParser;
p.StructExpand=true; % if we allow parameter  structure expanding
p.addRequired('I_org', @(x) isa(x,'uint8'));
p.addParamValue('tiledim', 30, @isnumeric); %Background correction. Overlap of 15 px, MinPts = 6 and eps = 0.1 are determined automatically
p.addParamValue('lambda', 5, @isnumeric); %Segmentation
p.addParamValue('minSizeMSER', 30, @isnumeric); %Segmentation
p.addParamValue('maxSizeMSER', 4000, @isnumeric); %Segmentation
p.addParamValue('maxVariation', 1, @isnumeric); %Segmentation
p.addParamValue('maxEcc', .7, @isnumeric); %Cell Splitting
p.addParamValue('minSizeSplit', 30, @isnumeric); %Cell Splitting
p.addParamValue('maxSizeSplit', 1000, @isnumeric); %Cell Splitting
p.addParamValue('visualize', false, @islogical); %Visualization

p.parse(I_org,varargin{:});
r = p.Results;

I_org = im2double(I_org);
%% compute the background
bg = bgest(I_org,r.tiledim);

%% correct the image
I = I_org./bg;
I(I>1) = 1;
I(I<0) = 0;

I=im2uint8(I);

%% segment the image
msers = linearMser(imcomplement(I),r.lambda,r.minSizeMSER,r.maxSizeMSER,r.maxVariation,0);
bw = zeros(size(I));
for mser=1:numel(msers)
    bw(msers{mser}) = 1;
end
bw = logical(bw);

%% split clumped cells
[L,bw2] = splitCells(I,bw,r.minSizeSplit,r.maxSizeSplit,r.maxEcc,1,1);

%% visualize it
if r.visualize
    s1 = subplot(2,2,1);
    imagesc(I_org)
    colormap gray
    title('Raw Image')
    axis off
    
    s2 = subplot(2,2,2);
    imagesc(I)
    colormap gray
    title('Background corrected')
    axis off
    
    s3 = subplot(2,2,3);
    imagesc(bw)
    colormap gray
    title('Segmented')
    axis off
    
    s4 = subplot(2,2,4);
    imagesc(bw2)
    colormap gray
    title('Final result after split/merge')
    axis off
    
    linkaxes([s1 s2 s3 s4])
end
end




% splitCells.m
function [L,bw2] = splitCells(I_c,bw,MinDiameter,MaxDiameter,maxEcc,docComplement,doMerge)

if docComplement
    I_c = imcomplement(I_c);
    I_c = im2double(I_c);
end

clearborder = 1;

bw = imfill(bw,'holes');

MaximaSuppressionSize = round(MinDiameter/1.5);

MaximaImage = bwulterode(bw);

DistanceTransformedImage = bwdist(~bw);

Overlaid = imimposemin(-DistanceTransformedImage,MaximaImage);

WatershedBoundaries = watershed(Overlaid) > 0;
bw = bw.*WatershedBoundaries;

bw = bwlabel(bw);

bw_woMerge = bwlabel(bw > 0);
bw_woMerge(bw_woMerge > 0) = 1;

if ~doMerge
    bw2 = bw_woMerge;
    L = bw2;
    return
end

%%% Remove bw with no marker in them (this happens occasionally)
%%% This is a very fast way to get pixel indexes for the bw
tmp = regionprops(bw,'PixelIdxList');
for k = 1:length(tmp)
    %%% If there is no maxima in these pixels, exclude object
    if sum(MaximaImage(tmp(k).PixelIdxList)) == 0
        bw(tmp(k).PixelIdxList) = 0;
    end
end


%%% Label the bw
bw = bwlabel(bw);

%%% Merge small bw

if doMerge
    NumberOfbwBeforeMerge = max(bw(:));
    bw = Mergebw(bw,I_c,[MinDiameter MaxDiameter],maxEcc);
    NumberOfbwAfterMerge = max(bw(:));
    NumberOfMergedbw = NumberOfbwBeforeMerge-NumberOfbwAfterMerge;
end

%%% Remove bw along the border of the image (depends on user input)
tmp = bw;
if clearborder
    bw = imclearborder(bw);
end

%%% Relabel the bw
[bw,NumOfbw] = bwlabel(bw > 0);
L = logical(bw);
bw2 = L;
bw2(L ~= 0) = 1;
end


function bw = Mergebw(bw,I_c,Diameters,MaxEccentricity)

%%% Find the object that we should try to merge with other bw. The object
%%% numbers of these bw are stored in the variable 'MergeIndex'. The bw
%%% that we will try to merge are either the ones that fall below the specified
%%% MinDiameter threshold, or relatively small bw that are above the MaxEccentricity
%%% threshold. These latter bw are likely to be cells where two maxima have been
%%% found and the watershed transform has divided cells into two parts.
MinDiameter = Diameters(1);
MaxDiameter = Diameters(2);

warning('off', 'MATLAB:divideByZero'); %%% Matlab failing atan vs atan2 in regionprops line 672.
props = regionprops(bw,'EquivDiameter','PixelIdxList','Eccentricity');   % Get diameters of the bw
warning('on', 'MATLAB:divideByZero');
EquivDiameters = cat(1,props.EquivDiameter);
Eccentricities = cat(1,props.Eccentricity);
IndexEccentricity = intersect(find(Eccentricities > MaxEccentricity),find(EquivDiameters < (MinDiameter + (MaxDiameter - MinDiameter)/4)));
IndexDiameter = find(EquivDiameters < MinDiameter);
MergeIndex = unique([IndexDiameter;IndexEccentricity]);

% Try to merge until there are no bw left in the 'MergeIndex' list.
[sr,sc] = size(I_c);
counter = round(numel(MergeIndex)/10);
prevcounter = counter;
while ~isempty(MergeIndex)
    counter = round(numel(MergeIndex)/10);
    if counter < prevcounter
        fprintf('MergeIndex: %i\n',numel(MergeIndex))
        prevcounter = counter;
    end
    
    % Get next object to merge
    CurrentObjectNbr = MergeIndex(1);
    
    %%% Identify neighbors and put their label numbers in a list 'NeighborsNbr'
    %%% Cut a patch so we don't have to work with the entire image
    [r,c] = ind2sub([sr sc],props(CurrentObjectNbr).PixelIdxList);
    rmax = min(sr,max(r) + 3);
    rmin = max(1,min(r) - 3);
    cmax = min(sc,max(c) + 3);
    cmin = max(1,min(c) - 3);
    bwPatch = bw(rmin:rmax,cmin:cmax);
    BinaryPatch = double(bwPatch == CurrentObjectNbr);
    GrownBinaryPatch = conv2(BinaryPatch,double(getnhood(strel('disk',2))),'same') > 0;
    Neighbors = bwPatch .*GrownBinaryPatch;
    NeighborsNbr = setdiff(unique(Neighbors(:)),[0 CurrentObjectNbr]);
    
    
    %%% For each neighbor, calculate a set of criteria based on which we decide if to merge.
    %%% Currently, two criteria are used. The first is a Likelihood ratio that indicates whether
    %%% the interface pixels between the object to merge and its neighbor belong to a background
    %%% class or to an object class. The background class and object class are modeled as Gaussian
    %%% distributions with mean and variance estimated from the image. The Likelihood ratio determines
    %%% to which of the distributions the interface voxels most likely belong to. The second criterion
    %%% is the eccentrity of the object resulting from a merge. The more circular, i.e., the lower the
    %%% eccentricity, the better.
    LikelihoodRatio    = zeros(length(NeighborsNbr),1);
    MergedEccentricity = zeros(length(NeighborsNbr),1);
    for j = 1:length(NeighborsNbr)
        
        %%% Get Neigbor number
        CurrentNeighborNbr = NeighborsNbr(j);
        
        %%% Cut patch which contains both original object and the current neighbor
        [r,c] = ind2sub([sr sc],[props(CurrentObjectNbr).PixelIdxList;props(CurrentNeighborNbr).PixelIdxList]);
        rmax = min(sr,max(r) + 3);
        rmin = max(1,min(r) - 3);
        cmax = min(sc,max(c) + 3);
        cmin = max(1,min(c) - 3);
        bwPatch = bw(rmin:rmax,cmin:cmax);
        I_cPatch = I_c(rmin:rmax,cmin:cmax);
        
        %%% Identify object interiors, background and interface voxels
        BinaryNeighborPatch      = double(bwPatch == CurrentNeighborNbr);
        BinaryObjectPatch        = double(bwPatch == CurrentObjectNbr);
        GrownBinaryNeighborPatch = conv2(BinaryNeighborPatch,ones(3),'same') > 0;
        GrownBinaryObjectPatch   = conv2(BinaryObjectPatch,ones(3),'same') > 0;
        Interface                = GrownBinaryNeighborPatch.*GrownBinaryObjectPatch;
        Background               = ((GrownBinaryNeighborPatch + GrownBinaryObjectPatch) > 0) - BinaryNeighborPatch - BinaryObjectPatch - Interface;
        WithinObjectIndex        = find(BinaryNeighborPatch + BinaryObjectPatch);
        InterfaceIndex           = find(Interface);
        BackgroundIndex          = find(Background);
        
        %%% Calculate likelihood of the interface belonging to the background or to an object.
        WithinObjectClassMean   = mean(I_cPatch(WithinObjectIndex));
        WithinObjectClassStd    = std(I_cPatch(WithinObjectIndex)) + sqrt(eps);
        BackgroundClassMean     = mean(I_cPatch(BackgroundIndex));
        BackgroundClassStd      = std(I_cPatch(BackgroundIndex)) + sqrt(eps);
        InterfaceMean           = mean(I_cPatch(InterfaceIndex)); %#ok Ignore MLint
        LogLikelihoodObject     = -log(WithinObjectClassStd^2) - (InterfaceMean - WithinObjectClassMean)^2/(2*WithinObjectClassStd^2);
        LogLikelihoodBackground = -log(BackgroundClassStd^2) - (InterfaceMean - BackgroundClassMean)^2/(2*BackgroundClassStd^2);
        LikelihoodRatio(j)      =  LogLikelihoodObject - LogLikelihoodBackground;
        
        %%% Calculate the eccentrity of the object obtained if we merge the current object
        %%% with the current neighbor.
        MergedObject =  double((BinaryNeighborPatch + BinaryObjectPatch + Interface) > 0);
        tmp = regionprops(MergedObject,'Eccentricity');
        MergedEccentricity(j) = tmp(1).Eccentricity;
        
        %%% Get indexes for the interface pixels in original image.
        %%% These indexes are required if we need to merge the object with
        %%% the current neighbor.
        tmp = zeros(size(I_c));
        tmp(rmin:rmax,cmin:cmax) = Interface;
        tmp = regionprops(double(tmp),'PixelIdxList');
        OrigInterfaceIndex{j} = cat(1,tmp.PixelIdxList); %#ok Ignore MLint
    end
    
    %%% Let each feature rank which neighbor to merge with. Then calculate
    %%% a score for each neighbor. If the neighbors is ranked 1st, it will get
    %%% 1 point; 2nd, it will get 2 points; and so on. The lower score the better.
    [ignore,LikelihoodRank]   = sort(LikelihoodRatio,'descend'); %#ok Ignore MLint % The higher the LikelihoodRatio the better
    [ignore,EccentricityRank] = sort(MergedEccentricity,'ascend'); %#ok Ignore MLint % The lower the eccentricity the better
    NeighborScore = zeros(length(NeighborsNbr),1);
    for j = 1:length(NeighborsNbr)
        NeighborScore(j) = find(LikelihoodRank == j) +  find(EccentricityRank == j);
    end
    
    %%% Go through the neighbors, starting with the highest ranked, and merge
    %%% with the first neighbor for which certain basic criteria are fulfilled.
    %%% If no neighbor fulfil the basic criteria, there will be no merge.
    [ignore,TotalRank] = sort(NeighborScore); %#ok Ignore MLint
    for j = 1:length(NeighborsNbr)
        CurrentNeighborNbr = NeighborsNbr(TotalRank(j));
        
        %%% To merge, the interface between bw must be more likely to belong to the object class
        %%% than the background class. The eccentricity of the merged object must also be lower than
        %%% for the original object.
        if LikelihoodRatio(TotalRank(j)) > 0 && MergedEccentricity(TotalRank(j)) < MaxEccentricity%Eccentricities(CurrentObjectNbr)
            
            %%% OK, let's merge!
            %%% Assign the neighbor number to the current object
            bw(props(CurrentObjectNbr).PixelIdxList) = CurrentNeighborNbr;
            
            %%% Assign the neighbor number to the interface pixels between the current object and the neigbor
            bw(OrigInterfaceIndex{TotalRank(j)}) = CurrentNeighborNbr;
            
            %%% Add the pixel indexes to the neigbor index list
            props(CurrentNeighborNbr).PixelIdxList = cat(1,...
                props(CurrentNeighborNbr).PixelIdxList,...
                props(CurrentObjectNbr).PixelIdxList,...
                OrigInterfaceIndex{TotalRank(j)});
            
            %%% Remove the neighbor from the list of bw to be merged (if it's there).
            MergeIndex = setdiff(MergeIndex,CurrentNeighborNbr);
        end
    end
    
    %%% OK, we are done with the current object, let's go to the next
    MergeIndex = MergeIndex(2:end);
end

%%% Finally, relabel the bw
bw = bwlabel(bw > 0);
end

function [xpts,ypts] = getpoints(AxisHandle)

Position = get(AxisHandle,'Position');
FigureHandle = (get(AxisHandle, 'Parent'));
PointHandles = [];
xpts = [];
ypts = [];
NbrOfPoints = 0;
done = 0;
%%% Turns off the CPimagetool function because it interferes with getting
%%% points.
ImageHandle = get(AxisHandle,'children');
set(ImageHandle,'ButtonDownFcn','');

hold on
while ~done;
    
    UserInput = waitforbuttonpress;                            % Wait for user input
    SelectionType = get(FigureHandle,'SelectionType');         % Get information about the last button press
    CharacterType = get(FigureHandle,'CurrentCharacter');      % Get information about the character entered
    
    % Left mouse button was pressed, add a point
    if UserInput == 0 && strcmp(SelectionType,'normal')
        
        % Get the new point and store it
        CurrentPoint  = get(AxisHandle, 'CurrentPoint');
        xpts = [xpts CurrentPoint(2,1)];
        ypts = [ypts CurrentPoint(2,2)];
        NbrOfPoints = NbrOfPoints + 1;
        
        % Plot the new point
        h = plot(CurrentPoint(2,1),CurrentPoint(2,2),'r.');
        set(AxisHandle,'Position',Position)                   % For some reason, Matlab moves the Title text when the first point is plotted, which in turn resizes the image slightly. This line restores the original size of the image
        PointHandles = [PointHandles h];
        
        % If there are any points, and the right mousebutton or the backspace key was pressed, remove a points
    elseif NbrOfPoints > 0 && ((UserInput == 0 && strcmp(SelectionType,'alt')) || (UserInput == 1 && CharacterType == char(8)))   % The ASCII code for backspace is 8
        
        NbrOfPoints = NbrOfPoints - 1;
        xpts = xpts(1:end-1);
        ypts = ypts(1:end-1);
        delete(PointHandles(end));
        PointHandles = PointHandles(1:end-1);
        
   elseif NbrOfPoints >= 3 && UserInput == 1 && CharacterType == char(13)
        done = 1;
        if ~isempty(PointHandles)
            delete(PointHandles)
        end
    end
end
xpts = round(xpts);
ypts = round(ypts);
hold off
set(ImageHandle,'ButtonDownFcn','CPimagetool');
end