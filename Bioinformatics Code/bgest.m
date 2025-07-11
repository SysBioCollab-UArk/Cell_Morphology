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