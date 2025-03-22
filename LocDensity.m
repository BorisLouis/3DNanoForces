%% Localization density


%% input
InputDirectory = 'F:\04-Apr\3 - Colloidal memory\3D DEP trapping\0.5V 100 kHz';

load([InputDirectory '\' 'trackResults.mat']);

TrackedData = trackRes.traces;

%% Filter out traces that are not going to the trap
for i=1:length(TrackedData)
    
    currTrace = TrackedData{i,1};
    %test= sum(gradient(currTrace.z));
%     test1 = sum(gradient(currTrace.z))>2000;
%     test2 = norm([currTrace.col(1), currTrace.row(1)]-CM)> norm([currTrace.col(end), currTrace.row(end)]-CM);
%     test3 = height(currTrace)>5;
%     test4 = and(currTrace.col(end)< CM(1)+2000,currTrace.col(end)> CM(1)-2000); 
%     test5 = and(currTrace.row(end)< CM(2)+2000,currTrace.row(end)> CM(2)-2000); 
%     test = logical(test1*test2*test3*test4*test5);
    
    if test
        figure(1)
        hold on
        
        plot3(currTrace.col,currTrace.row,currTrace.z);
        axis image
    else
        TrackedData{i,1} = [];
    end
    
    
end
%%
traces2Keep = TrackedData;
idx = cellfun(@isempty,traces2Keep);
traces2Keep(idx,:) = [];
%% Extract all cooredinate
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
   % traces2Keep{i,1} = currTrace;
    
    allRow = [allRow ;currTrace.row];
    allCol = [allCol ; currTrace.col];
    allz = [allz ; currTrace.z];
    allu = [allu ; currTrace.u];
    allv = [allv; currTrace.v];
    allw = [allw; currTrace.w];

    allSpeed = [allSpeed ;currTrace.speedR];
    maxSpeed(i) = max(currTrace.speedR);
    
end

%% binning
step = 200;
ybin = 10000:step:17000;
xbin = 18000:step:25000;
zStep =-200;
counter = 1;
maxCount = 50;
zVec = 2500:zStep:-1800;
figure
hold on
im = zeros(length(ybin)-1,length(xbin)-1,length(zVec)-1);
BW = im;

for z = zVec
    currentCol = allCol(allz<z & allz>z+zStep);
    currentRow = allRow(allz<z & allz>z+zStep);
    [N,Xedges,Yedges] = histcounts2(currentCol,currentRow,xbin,ybin);
    im(:,:,counter) = N;
    
    tmpBW = N>0;
     tmpBW = bwareaopen(tmpBW,10);
      se = strel('disk',3);
      tmpBW = imclose(tmpBW,se);
%     
    BW(:,:,counter) = tmpBW;
    
    
%     subplot(3,length(zVec)/3,counter)
%     imagesc(Xedges,Yedges,N)
%     title(['z = ', num2str(z)]);
%     axis image
%     colormap('hot')
%     caxis([0 maxCount])
    counter= counter+1;
    
end
BW(:,:,end) = zeros(size(BW(:,:,end)));
% Displaying 3D model
mex +rendering3D\smoothpatch_curvature_double.c -v
mex +rendering3D\smoothpatch_inversedistance_double.c -v
mex +rendering3D\vertex_neighbours_double.c -v

data2Render = flip(BW,3);
iSurface = isosurface(data2Render,1/2);

% smoothing using compiled c code
smoothISurface = rendering3D.smoothpatch(iSurface,0,10);
%comnvert to px
smoothISurface.vertices(:,1) = (smoothISurface.vertices(:,1));
smoothISurface.vertices(:,2) = (smoothISurface.vertices(:,2));
smoothISurface.vertices(:,3) = (smoothISurface.vertices(:,3));


%z-coloring
colorModel = smoothISurface.vertices(:,3)/max(smoothISurface.vertices(:,3));
zColor = true;

%Plot the network with Z coloring or unique color depending on the user
%input
figure
if zColor
    p = patch('Faces',smoothISurface.faces,'Vertices',smoothISurface.vertices,'FaceVertexCData',colorModel,'FaceColor','interp');
    colormap('jet')
    p.EdgeColor = 'none';
    daspect([2 2 1])
    view(3);
    axis tight
    camlight
    lighting gouraud
    title('Z-coloring')
else
    p2 = patch(smoothISurface);
    p2.FaceColor = colorModel;
    p2.EdgeColor = 'none';
    view(3);
    axis tight
    camlight
    lighting gouraud
    title('unicolor');
end
%% binarizing image




