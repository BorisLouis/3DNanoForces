%store all the data from different movies in one large array
nData =sum(cellfun(@height,trackRes.traces));
allData = zeros(nData(1),4);
for i = 1:size(trackRes.traces,1)
    coin = rand(1);
    currTrace = trackRes.traces{i,1};
    
    id = find(allData(:,1)==0,1,'first');
    allData(id:id+height(currTrace)-1,1) = currTrace.col;
    allData(id:id+height(currTrace)-1,2) = currTrace.row;
    allData(id:id+height(currTrace)-1,3) =currTrace.z;
    allData(id:id+height(currTrace)-1,4) =currTrace.rT/1000;
       
end

%% Calculate 3D Histogram
%step size for histogram sampling in 3D
stepZ = 800;
stepXY = 200;
zRange = max(allData(:,3)) - min(allData(:,3));
xRange = max(allData(:,1)) - min(allData(:,1));
yRange = max(allData(:,2)) - min(allData(:,2));
%make the 3D grid to bin the data
edgesZ = linspace(min(allData(:,3)),max(allData(:,3)),round(zRange/stepZ));
edgesX = linspace(0,max(allData(:,1)),round(xRange/stepXY));
edgesY = linspace(0,max(allData(:,2)),round(yRange/stepXY));
%keep track of the bins
binZ = edgesZ(2)-edgesZ(1);
binY = edgesY(2)-edgesY(1);
binX = edgesX(2)-edgesX(1);
% edgesZ = min(allData(:,3)):stepZ:max(allData(:,3));
% edgesX = min(allData(:,1)):stepXY:max(allData(:,1));
% 
% edgesY = min(allData(:,2)):stepXY:max(allData(:,2));

%binning, for speed we use matlab binning in2D after preselecting the data
%based on Z, diminishing greatly the number of iteration
histogram3D = zeros(length(edgesY)-1,length(edgesX)-1,length(edgesZ)-1);

for i = 1:length(edgesZ)-1
    
    currData = allData(and(allData(:,3)>=edgesZ(i), allData(:,3) <edgesZ(i+1)),:);

    [histogram3D(:,:,i)] = histcounts2(currData(:,2),currData(:,1),edgesY,edgesX);
    
end
%%
%plot the last frame (interface)
figure
imagesc(histogram3D(:,:,end))
axis image
clim([0 30])
title([num2str(edgesZ(end-1)) '-' num2str(edgesZ(end)) '\nm'])

hold on

%% Plot all Traces with time color-coding (4D plot)
minSize = 50;%number of frame the traces needs to last to be plotted.
expTime = 0.008; %sec
pxSize = 200;

maxFr = zeros(size(trackRes.traces,1),1);
for i = 1:size(trackRes.traces,1)
    
    id = round(rand(1)*size(trackRes.traces,1));

    currTrace = trackRes.traces{id,1};
    
    if height(currTrace) > minSize
        colPlot = currTrace.col/pxSize;
        rowPlot = currTrace.row/pxSize;
        zPlot   = currTrace.z;
        tPlot   = currTrace.t*expTime;
        pl = plot(colPlot,rowPlot,'r');
        %plot with time color coding
        %patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')
        disp(id);  
        w = waitforbuttonpress;
         delete(pl);
         %disp(id);   
    end
    maxFr(i,:) = max(currTrace.t);
    
end

maxFr = max(maxFr);
view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image


%% plot specific trace

id = 9122;

%%
figure (1)
clf
id = 19969;
subplot(1,2,1)
imagesc(histogram3D(:,:,end))
hold on
currTrace = trackRes.traces{id,1};
colPlot = currTrace.col/200;
rowPlot = currTrace.row/200;
plot(colPlot,rowPlot,'r');

subplot(1,2,2)
currTrace = trackRes.traces{id,1};
colPlot = currTrace.col;
rowPlot = currTrace.row;
zPlot   = currTrace.z;
tPlot   = (currTrace.t-currTrace.t(1))*expTime;
%pl = plot(colPlot,rowPlot,'r');
%plot with time color coding
patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')
axis image
box on 
colorbar
view(3)

figure(2)
clf
patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')
axis image
box on 
colorbar
view(3)



% for i = 1:size(trackRes.traces,1)
% 
%     currTrace = trackRes.traces{i,1};
% 
%     if and(height(currTrace)>20,any(and(24000<round(currTrace.row), round(currTrace.row)< 24400)))
% 
%         disp(i);
% 
%     end
% 
% end


%% plot traces on top of electrodes for SI
%id1 ejection perpendicular
list.id1 = [83,8550, 2886,17014,15429,3650,11446, 23885, 6227, 6909,19969];
%id2 Trapping event
list.id2 = [26606,967, 31037, 14425, 24326,31246, 31073];
% id 3 trapped
list.id3 = [9688,11764,3340,11267];
% id 4 vortex
list.id4 = [11790, 22763,15870, 6160, 22763, 25952,9323]; 
%id5 Brownian
list.id5 = [11801, 16002];
%id6 = 45 degrees
list.id6 = [32210, 21668,30453,22930,30577,2244,2496,7167,9683,16738 ];

col = [0.3010 0.7450 0.9330;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.6350 0.0780 0.1840];

figure
hold on
imagesc(histogram3D(:,:,end))
colormap('hot')
for i=1:6
    currentList = ['id' num2str(i)];

    data = list.(currentList);

    for j = 1:length(data)
        currTrace = trackRes.traces{data(j),1};
        colPlot = currTrace.col/200;
        rowPlot = currTrace.row/200;
        plot(colPlot,rowPlot,'color',col(i,:));

    end
end
axis image
xlim([1 159])
ylim([1 250])

