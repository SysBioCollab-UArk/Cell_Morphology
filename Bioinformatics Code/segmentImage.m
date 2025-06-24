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