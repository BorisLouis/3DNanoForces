%% input
InputDirectory = 'D:\Documents\Unif\PhD\2020-Data\01-Jan\14 - Trapping\TrappingData\1000mW - 2';

load([InputDirectory '\' 'trackResults.mat']);

TrackedData = trackRes.traces;
CM = [23284,10027];
expTime = 0.005;

color = getColorFromCmap('lines',length(TrackedData));

%%
for i=1:length(TrackedData)
    
    currTrace = TrackedData{i,1};
    %test= sum(gradient(currTrace.z));
    test1 = sum(gradient(currTrace.z))>2000;
    test2 = norm([currTrace.col(1), currTrace.row(1)]-CM)> norm([currTrace.col(end), currTrace.row(end)]-CM);
    test3 = height(currTrace)>5;
    test4 = and(currTrace.col(end)< CM(1)+2000,currTrace.col(end)> CM(1)-2000); 
    test5 = and(currTrace.row(end)< CM(2)+2000,currTrace.row(end)> CM(2)-2000); 
    test = logical(test1*test2*test3*test4*test5);
    
    if test
        figure(1)
        hold on
        
        plot3(currTrace.col,currTrace.row,currTrace.z,'color',color(i,:));
        axis image
    else
        TrackedData{i,1} = [];
    end
    
    
end
%%
traces2Keep = TrackedData(~cellfun(@isempty,TrackedData(:,1)),1);
allRow =[];
allSpeed =[];
allz =[];
allCol = [];
allu = [];
allv = [];
allw = [];
maxSpeed = zeros(1,length(traces2Keep));

for i=1:length(traces2Keep)
    currTrace = traces2Keep{i,1};
    
    for j=2:height(currTrace)-1
        
        currTrace.u(j) = (diff(currTrace.col([j-1 j+1])));
        currTrace.v(j) = (diff(currTrace.row([j-1 j+1])));
        currTrace.w(j)   = (diff(currTrace.z([j-1 j+1])));
        currTrace.speedR(j)   = sqrt(currTrace.u(j).^2+currTrace.v(j).^2+currTrace.w(j).^2);
       

    end

    currTrace.speedR(currTrace.speedR==0) = NaN;
    currTrace.u(currTrace.u==0) = NaN;
    currTrace.v(currTrace.v==0) = NaN;
    currTrace.w(currTrace.w==0) = NaN;
    traces2Keep{i,1} = currTrace;
    
    allRow = [allRow ;currTrace.row];
    allCol = [allCol ; currTrace.col];
    allz = [allz ; currTrace.z];
    allu = [allu ; currTrace.u];
    allv = [allv; currTrace.v];
    allw = [allw; currTrace.w];

    allSpeed = [allSpeed ;currTrace.speedR];
    maxSpeed(i) = max(currTrace.speedR);
    
end

disp(['average speed: ' num2str(nanmax(allSpeed))])

%% Plot Image

figure(2)
hold on
quiver3(allCol,allRow,allz,allu,allv,allw);


%%

currRow = allRow(allz<-3500);
currCol = allCol(allz<-3500);

%[N,rowEdges,zEdges] = histcounts2(currRow,allz,500);

[N,rowEdges,colEdges] = histcounts2(currRow,currCol,200);


speedBinned = rowEdges;
for i=1:length(rowEdges)-1
    for j =1:length(colEdges)-1
        
        speedBinned(i,j) = nanmean(allSpeed(currRow<rowEdges(i+1) & currRow>=rowEdges(i)...
            & currCol<colEdges(j+1) & currCol>=colEdges(j)));
    
    end
end

figure
imagesc(colEdges,rowEdges,speedBinned)
caxis([0 200])
axis image
