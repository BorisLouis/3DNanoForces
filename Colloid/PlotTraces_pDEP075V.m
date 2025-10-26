

% path2Save = 'H:\11032025 vortex testing\L3_4W\10ms';
% ext = '.gif';
% filename=sprintf('%s%sdata%s', path2Save,filesep,ext);

minSize = 3;%number of frame the traces needs to last to be plotted.
expTime = 0.01; %sec
sizeParticles = 500; % diameter in nm
frameRate = 100;
trailing = 20; %frame the traces stays in the movie

perc = 0.4;%proportion of traces to plot

%% Top View with time color-coding (4D plot)
CM = zeros(size(trackRes.traces,1),3);
maxFr = zeros(size(trackRes.traces,1),1);
figure
hold on
for i = 1:size(trackRes.traces,1)
    coin = rand(1);
    currTrace = trackRes.traces{i,1};
    
    if and(height(currTrace) > minSize,coin<=perc)
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        tPlot   = (currTrace.t-currTrace.t(1))*expTime;

        tPlot   = tPlot/(max(tPlot));
        plot3(colPlot,rowPlot,zPlot)
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

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
ylim(xlim+10^4)
view(2)
box on
set(gcf,'color','w')

%% %% Top View with time color-coding (4D plot)
CM = [24600 15770];
perc = 0.4;
figure
hold on
for i = 1:size(trackRes.traces,1)
    coin = rand(1);
    currTrace = trackRes.traces{i,1};
    
    if and(height(currTrace) > minSize,coin<=perc)
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        tPlot   = (currTrace.t-currTrace.t(1))*expTime;

        tPlot   = tPlot/(max(tPlot));
        plot3(colPlot,rowPlot,zPlot)
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

    end        
    % CM(i,:) = [mean(currTrace.row),mean(currTrace.col),mean(currTrace.z)];
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
r = 5000;
xlim([CM(2)-r,CM(2)+r])
ylim([CM(1)-r,CM(1)+r])
view(3)
box on
set(gcf,'color','w')

%% %% Cross View with time color-coding (4D plot)

perc = 0.4;
figure
hold on
for i = 1:size(trackRes.traces,1)
    coin = rand(1);
    currTrace = trackRes.traces{i,1};
    
    if and(height(currTrace) > minSize,coin<=perc)
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        tPlot   = (currTrace.t-currTrace.t(1))*expTime;

        tPlot   = tPlot/(max(tPlot));
        plot3(colPlot,rowPlot,zPlot)
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

    end        
    % CM(i,:) = [mean(currTrace.row),mean(currTrace.col),mean(currTrace.z)];
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
r = 3000;
xlim([CM(2)-5000-r,CM(2)-5000+r])
ylim([CM(1)-5000-r,CM(1)-5000+r])
view(135.8,-0)
%view(gca,[90 0])
box on
set(gcf,'color','w')