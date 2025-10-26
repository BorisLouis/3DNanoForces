%% 
% This code assumes that some trackRes has been loaded onto matlab
% no clearing
clc
%close all

%% User input

path2Save = 'P:\20250620_notrap_ps_motion\Conc._diffusion\100LYS\Pending';
ext = '.gif';
filename=sprintf('%s%sdata%s', path2Save,filesep,ext);

minSize = 5;%number of frame the traces needs to last to be plotted.
expTime = 0.05; %sec
sizeParticles = 200; % diameter in nm
frameRate = 100;
trailing = 20; %frame the traces stays in the movie


%% get Mask based on localization density





%% Plot all Traces with time color-coding (4D plot)
CM = zeros(size(trackRes.traces,1),3);
maxFr = zeros(size(trackRes.traces,1),1);
figure
hold on
for i = 1:size(trackRes.traces,1)
    
    currTrace = trackRes.traces{i,1};
    
    if height(currTrace) > minSize
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        t = (currTrace.t-currTrace.t(1));
        tPlot   =t ./max(t);
        zT = (currTrace.z-currTrace.z(1));
   
        zCol = zT ./max(zT);
        
        if and(any(zPlot>1800),any(diff(colPlot)<50))

            %plot3(colPlot,rowPlot,zPlot)
            %plot with time color coding
            patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')
        end
    end        
    CM(i,:) = [mean(currTrace.row),mean(currTrace.col),mean(currTrace.z)];
    maxFr(i,:) = max(currTrace.t);
end
CM = mean(CM,1);
maxFr = max(maxFr);
view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
set(gcf,'color','w')


%% test Plot all Traces with time color-coding (4D plot)
CM = zeros(size(trackRes.traces,1),3);
maxFr = zeros(size(trackRes.traces,1),1);
figure
hold on
i = 2
    
    currTrace = trackRes.traces{i,1};
    
    if height(currTrace) > minSize
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        tPlot   = currTrace.t*expTime;
        plot3(colPlot,rowPlot,zPlot)
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

    end        
    CM(i,:) = [mean(currTrace.row),mean(currTrace.col),mean(currTrace.z)];
    maxFr(i,:) = max(currTrace.t);

CM = mean(CM,1);
maxFr = max(maxFr);
view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
%% Make Awesome movie

radius = 1000;
yLimit = [0 35000];
xLimit = [0 25000];
zLimit = [-2000 2000];


Fig = figure;
view(3);
xlim(xLimit);
ylim(yLimit);
zlim(zLimit);
xlim manual;
ylim manual;
gcf;
hold on

[x,y,z] = sphere(32);
x = x*sizeParticles/2;
y = y*sizeParticles/2;
z = z*sizeParticles/2;

camlight
lighting('gouraud');
for i = 1 :maxFr
    for j = 1:size(trackRes.traces,1)
        currTrace = trackRes.traces{j,1};
        if height(currTrace) > minSize
            idx2Frame = currTrace.t==i;
            idx = i-trailing:i;
            idx = ismember(currTrace.t,idx);
            if and(~all(idx==0), ~all(idx2Frame ==0))

                data2Plot = currTrace(idx,:);
                axis image
                xlim(xLimit);
                ylim(yLimit);
                zlim(zLimit)
                xlim manual;
                ylim manual;
                zlim manual;

                view(3);

                gcf;
                hold on

                plot3(data2Plot.col,data2Plot.row,data2Plot.z,'color',[0 0 1])

                X = x+currTrace.col(idx2Frame);
                Y = y+currTrace.row(idx2Frame);
                Z = z+currTrace.z(idx2Frame);

                surf(X,Y,Z,'LineStyle','none','Facecolor',[0.4,0.4,0.4])
                

            end
        end
    end
    drawnow;
    frame = getframe(Fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if i == 1

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

    else

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

    end
    clf;

end

%% plot only a single traced
figure
idx = 10;

currTrace = trackRes.traces{idx,1};
    
colPlot = currTrace.col;
rowPlot = currTrace.row;
zPlot   = currTrace.z;
tPlot   = currTrace.t*expTime;

patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image