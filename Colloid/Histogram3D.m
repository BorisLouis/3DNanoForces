%%
cleanUpStillParticles =false;

%% 3D Histogram
nData =sum(cellfun(@height,trackRes.traces));
allData = zeros(nData(1),3);
for i = 1:size(trackRes.traces,1)
    coin = rand(1);
    currTrace = trackRes.traces{i,1};
    
    id = find(allData(:,1)==0,1,'first');
    allData(id:id+height(currTrace)-1,1) = currTrace.col;
    allData(id:id+height(currTrace)-1,2) = currTrace.row;
    allData(id:id+height(currTrace)-1,3) =currTrace.z;      
       
end

%% Histogram 3d
stepZ = 1000;
stepXY = 200;

% edgesZ = linspace(min(allData(:,3)),max(allData(:,3)),stepZ);
% edgesX = linspace(min(allData(:,1)),max(allData(:,2)),stepXY);
% edgesY = linspace(min(allData(:,2)),max(allData(:,3)),stepXY);


edgesZ = min(allData(:,3)):stepZ:max(allData(:,3));
edgesX = min(allData(:,1)):stepXY:max(allData(:,1));

edgesY = min(allData(:,2)):stepXY:max(allData(:,2));

histogram3D = zeros(length(edgesY)-1,length(edgesX)-1,length(edgesZ)-1);

for i = 1:length(edgesZ)-1
    
    currData = allData(and(allData(:,3)>=edgesZ(i), allData(:,3) <edgesZ(i+1)),:);

    [histogram3D(:,:,i)] = histcounts2(currData(:,2),currData(:,1),edgesY,edgesX);
    
end

%% 


figure
for i = 1:size(histogram3D,3)
    
    imagesc(histogram3D(:,:,i))
    axis image
    clim([0 50])
    title([num2str(edgesZ(i)) '-' num2str(edgesZ(i+1)) '\nm'])
    w = waitforbuttonpress;
    clf

end

%%
if cleanUpStillParticles
%histogram3D(histogram3D(:)>200) =0;
%detect and cleanup standing particles
bw = imbinarize(histogram3D(60:417,:,end),250);

bw = imdilate(bw,strel('disk',2));
else
    bw = zeros(size(histogram3D(60:417,:,end)));
end
figure
imagesc(histogram3D(50:407,:,end).*~bw)
axis image
clim([0 10])
colormap('parula')
colorbar


%% rotate histogram for averaging
if cleanUpStillParticles
rotHist = imrotate(histogram3D(60:417,:,:).*~bw,44.85);
else
    rotHist = imrotate(histogram3D(60:417,:,:),44.85);
end
figure
imagesc(rotHist(:,:,end))
axis image
clim([0 40])
colormap('parula')

figure
subplot(1,2,1)
imagesc(squeeze(sum(rotHist,1))')
clim([0 40])
subplot(1,2,2)
imagesc(squeeze(sum(rotHist,2))')
clim([0 100])



%%

CM = mean(CM,1);
maxFr = max(maxFr);
view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
ylim(xlim+10^4)
view(2)
box on
set(gcf,'color','w')