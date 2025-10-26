%% User Input
clear; clc;
Data = 'F:\WoutVT\20251003_trapping_statisics_BE_only_2 (1 MHz 2V)\200nm_PS_1MHz_2V';
folderControl = 'F:\WoutVT\20251003_trapping_statisics_BE_only_2 (1 MHz 2V)\200nm_PS_1MHz_Diffusion';

%Move parameters here

expTime = 0.01; %exposure time in seconds
zcond=1700;
minpoints=10;
%% get folder
C0 = dir(folderControl);

C0 = C0(3:end);

C0= C0([C0.isdir]);

%% Initial concentration
% Here we calculate concentration of particle, for this we need the
% dimension of the image:
vol = [383 556 8];

pixSize = 0.095;%size of pixels in XY
deltaZ = 0.580; % distance between planes in Z

volume = prod([vol(1)*pixSize,vol(2)*pixSize, vol(3)*deltaZ]);

% get concentration
% look at the number of particles detected in the volume
nPart = zeros(length(C0),1);
stdPart =nPart;
for i = 1:length(C0)
    currentFile = [C0(i).folder filesep C0(i).name filesep 'particle.mat'];
    
    currentPart = load(currentFile);
    %currentPart.particle.nParticles holds the number of particle for the
    %specific frame so we take the average of that number ot have the
    %average particle in the movie, we divide by volum and 1000 to get the
    % per 1000 micrometer cubes unit
    nPart(i) = mean((currentPart.particle.nParticles)/(volume/1000));
    stdPart(i) = std((currentPart.particle.nParticles)/(volume/1000));

end

figure

errorbar(1:length(C0),nPart,stdPart,'o','MarkerFaceColor','auto')
title('Concentration per 1000 \mum^3')
xlabel('Experiment cycle')
ylabel('Concentration in particle per 1000 \mum^3')


%% Localization density and concentration
%load the tracking file
load([Data filesep 'trackResults.mat'])

%% Conversion to array
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

%plotting the histogram
h = figure;
for i = 1:size(histogram3D,3)
    
    imagesc(histogram3D(:,:,i))
    axis image
    clim([0 30])
    title([num2str(edgesZ(i)) '-' num2str(edgesZ(i+1)) '\nm'])
    w = waitforbuttonpress;
    clf

end
close(h);

%plot the last frame (interface)
figure
 imagesc(histogram3D(:,:,end))
axis image
clim([0 30])
title([num2str(edgesZ(end-1)) '-' num2str(edgesZ(end)) '\nm'])


%% Detection of trapping event - Part 1
% Here we first do a segmentation of the area where there is significant
% number of localization as this should be the trapping area. %we transform
% it into a binary mask where all numbers are either 1 (there is a trapping
% site) or 0 (there is not trapping site) This will help us deciding when a
% particle was trapped or not).

Mask = histogram3D(:,:,end);

Mask = imbinarize(Mask,15);

Mask = bwareaopen(Mask,4);

se = strel('disk',3);

Mask = imdilate(Mask,se);

%% Detection of trapping event Part2 and speed calculations

trapEvent =table(zeros(1000,1),zeros(1000,1),'VariableNames',{'idx','time'});
figure
hold on
n= 1;
vecSpeed = []; 
FullPos = [];

%we loop through all traces to find which one ended up being trapped.
for i=1:length(trackRes.traces)
    currTrace = trackRes.traces{i};

    if (height(currTrace)>minpoints) %impose trace to have at least 4 datapoints
        
        currCol = currTrace.col;
        currRow = currTrace.row;
        %convert trace XY position to indices so we can compare it with our
        %mask of trapping site and know which traces landed in a trapping
        %site
        currInds = [floor((currCol)/binX)-1,floor((currRow)/binY)-1];
    
        inds = sub2ind(size(Mask),currInds(:,2),currInds(:,1));
    
        currZ   = currTrace.z;
        currT   = currTrace.rT/1000;
        %here comparison of indices with trapping site occurs
        test1 = Mask(inds);

        %Then we look if the z-position is high enough to be considered
        %trapped, current 1400 is considered high enough, can be changed
        
        if ~isempty(currZ)
            test2 = currZ(end)>zcond;% we make sure particles are in Z above
            %1400 in their last frame (they stay trapped) 
            %we test here that particle were not filling the condition at
            %the beginning of the measurement (otherwise was already stuck)
            test3=(and(test1(1)==1,currZ(1)>zcond));
        else
            test2=0;
            test3=0;
        end
       
        %we can compare both condition to get the one that were part of it        
        test = and(test1,test2);
        
        %here we count the event s and calculate speed
        if and(~all(test==0),~test3) %check if the particle is indeed in the criterion
            if length(test)>minpoints %we want the condition to be fulfilled for at least 4 frames
                %this is to ensure that it is not just an error particle
                %was in trap condition for a few frames 

                % here we get the timing of the event (when the condition
                % were first fulfilled
                eventTime = currT(test);
                trapEvent.idx(n) = i;
                trapEvent.time(n) = eventTime(1);
                n=n+1;
                %this is not currently used
                tPlot   = (currT-currT(1))*expTime;
                tPlot   = tPlot/(max(tPlot));
                
                %simple speed calculation
                ax = diff(currCol)/1000/expTime;
                ay = diff(currRow)/1000/expTime;
                az = diff(currZ)/1000/expTime;
                %get the speed distribution
                vecSpeed = [ vecSpeed; sqrt(ax.^2+ay.^2+az.^2)];

                %shift coordinate on top of each otherthem on top of each other:
                currCol = currCol-currCol(end);
                currRow = currRow-currRow(end);
                
                plot3(currCol,currRow,currZ)
        
                %plot with time color coding
                patch([currCol(:)' nan],[currRow(:)' nan],[currZ(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')
            
                x = currCol;
                y = currRow;
                z = currZ;
                FullPos = [FullPos; [ax,ay,az , x(1:end-1)+0.5, y(1:end-1)+0.5,z(1:end-1)+0.5]];
    
            
            end
        end

    end
end

view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image

view(2)
box on
set(gcf,'color','w')

%plot the trapping events
[N,Edges] = histcounts(trapEvent.time,1:30);
figure
subplot(1,2,1)
scatter(Edges(1:end-1)-0.5,N/10,20,"filled");
xlabel('Time(s)')
ylabel('Trap event')
axis square
box on
%title('Events as a function of time')
subplot(1,2,2)
scatter(Edges(1:end-1)-0.5,cumsum(N/10),20,"filled");
xlabel('Time(s)')
ylabel('Cumulative Trap event')
%title('Cumulative distribution of event as function of time')
axis square
box on

%plot the speed distribution
figure
histogram(vecSpeed,100)

%% vector plot
threshold =50;
filterPos =FullPos; %sum of XYZ should be >100
filterPos(sum(abs(filterPos(:,1:3)),2)<threshold,:)=[];

close all

cst = 20;
figure
quiverC3D(filterPos(:,4),filterPos(:,5),filterPos(:,6),...
    filterPos(:,1)*cst,filterPos(:,2)*cst,filterPos(:,3)*cst)
axis image
view(3)
colorbar
colormap('jet')
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image

%% binning for 2D vector plot ZX
%==> this would be better to do in spherical coordinate
binX = linspace(round(min(filterPos(:,4))),round(max(filterPos(:,4))),40);
binZ = linspace(round(min(filterPos(:,6))),round(max(filterPos(:,6))),40);

[X,Z] = meshgrid(binX(1:end-1),binZ(1:end-1));

U =zeros(size(X));
V = zeros(size(Z));
for i =1:length(binX)-1
    for j =1:length(binZ)-1

        cIdx = and(filterPos(:,4)>binX(i), filterPos(:,4)<=binX(i+1));
        cIdz  = and(filterPos(:,6)>binZ(j), filterPos(:,6)<=binZ(j+1));
        cId = and(cIdx,cIdz);
        U(j,i) = mean(nonzeros(filterPos(cId,1)));
        V(j,i) = mean(nonzeros(filterPos(cId,3)));
        
    end
end
% 
% figure
% quiver(X,Z,U,V)
cst = 15;

figure
quiverC3D(X,Z,zeros(size(X)),U*cst,V*cst,zeros(size(Z))*cst)
axis image
view(2)
colorbar
colormap('jet')

xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
xlim([-2000 2000])

figure
subplot(1,2,1)
streamline(X,Z,U,V,X,Z)
subplot(1,2,2)
streamslice(X,Z,U,V)
%% %% binning for 2D vector plot XY
%==> this would be better to do in spherical coordinate
binX = linspace(round(min(filterPos(:,4))),round(max(filterPos(:,4))),40);
binY = linspace(round(min(filterPos(:,5))),round(max(filterPos(:,5))),40);

[X,Y] = meshgrid(binX(1:end-1),binY(1:end-1));

U =zeros(size(X));
V = zeros(size(Y));
for i =1:length(binX)-1
    for j =1:length(binY)-1

        cIdx = and(filterPos(:,4)>binX(i), filterPos(:,4)<=binX(i+1));
        cIdz  = and(filterPos(:,5)>binY(j), filterPos(:,5)<=binY(j+1));
        cId = and(cIdx,cIdz);
        U(j,i) = mean(nonzeros(filterPos(cId,1)));
        V(j,i) = mean(nonzeros(filterPos(cId,3)));
        
    end
end
% 
% figure
% quiver(X,Z,U,V)
cst = 15;

figure
quiverC3D(X,Y,zeros(size(X)),U*cst,V*cst,zeros(size(Y))*cst)
axis image
view(2)
colorbar
colormap('jet')

xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
xlim([-2000 2000])
ylim([-2000 2000])
% figure
% subplot(1,2,1)
% streamline(X,Y,U,V,X,Y)
% subplot(1,2,2)
% streamslice(X,Y,U,Y)

%% Temporal concentration
%4D temporal concentration plot
stepZ = 400;
stepXY = 200;
stepT = 1; %100 = 1s at 10ms exposure time
% edgesZ = linspace(min(allData(:,3)),max(allData(:,3)),stepZ);
% edgesX = linspace(min(allData(:,1)),max(allData(:,2)),stepXY);
% edgesY = linspace(min(allData(:,2)),max(allData(:,3)),stepXY);

edgesT = min(allData(:,4)):stepT:max(allData(:,4));
edgesZ = min(allData(:,3)):stepZ:max(allData(:,3));
edgesX = min(allData(:,1)):stepXY:max(allData(:,1));
edgesY = min(allData(:,2)):stepXY:max(allData(:,2));

histogram4D = zeros(length(edgesY)-1,length(edgesX)-1,length(edgesZ)-1,length(edgesT));

for j=1:length(edgesT)-1
    for i = 1:length(edgesZ)-1
        
        currData = allData(and( ...
            and(allData(:,3)>=edgesZ(i), allData(:,4)>=edgesT(j)),...
            and(allData(:,3) <edgesZ(i+1),allData(:,4)<edgesT(j+1))),:);
    
        [histogram4D(:,:,i,j)] = histcounts2(currData(:,2),currData(:,1),edgesY,edgesX);
        
    end
end


figure
hold on
color = colormap("parula");

vec = round(linspace(1,length(color),size(histogram4D,4)-8));

col = color(vec,:);
for i = 1:size(histogram4D,4)-8
    
    data = squeeze(sum(sum(histogram4D(:,:,:,i),1),2))/100/10;%divided by the number of time unit per bin and the 10 movies

    plot(edgesZ(1:end-1),data,'-o','color',col(i,:))
    

end
colorbar
axis square
box on
xlabel('z position')
ylabel('Concentration')