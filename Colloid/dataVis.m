expTime = 0.008; %sec
minSize = 3; % minimum size for a trace

%% Plot all Traces with time color-coding (4D plot)

figure
hold on

for i = 1:size(trackRes.traces,1)
    
    currTrace = trackRes.traces{i,1};
    
    if height(currTrace) > minSize
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        %tPlot   = currTrace.t*expTime;
        tPlot   = (currTrace.t-currTrace.t(1))*expTime;
        tPlot   = tPlot/(max(tPlot));
        %plot3(colPlot,rowPlot,zPlot)
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

    end        
end

view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image