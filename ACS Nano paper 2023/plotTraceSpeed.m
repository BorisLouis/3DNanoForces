
currTrace = traces2Keep{149,1};
expTime = 0.005;

colPlot = currTrace.col;
rowPlot = currTrace.row;
zPlot   = currTrace.z;
tPlot   = currTrace.t*expTime;
tPlot = tPlot-tPlot(1);
%plot3(colPlot,rowPlot,zPlot)
%plot with time color coding
patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')
axis image
view(3)

colPlot = currTrace.col-CM(1);
rowPlot = currTrace.row-CM(2);
zPlot   = currTrace.z -CM(3);

rPlot = sqrt((colPlot).^2 + (rowPlot).^2 + (zPlot).^2);
speed = diff(rPlot)/1000/0.005;
speed = (diff(colPlot)+ diff(rowPlot) + diff(zPlot))/1000/0.005;


allSpeed = [];
CM = [23284,10027, 0];
for i = 1 : size(traces2Keep,1)
    
    currTrace = traces2Keep{i,1};
    expTime = 0.005;

    colPlot = currTrace.col-CM(1);
    rowPlot = currTrace.row-CM(2);
    zPlot   = currTrace.z -CM(3);
    
    rPlot = sqrt((colPlot).^2 + (rowPlot).^2 + (zPlot).^2);
    
    zPlot   = currTrace.z;
    tPlot   = currTrace.t*expTime;
    tPlot = tPlot-tPlot(1);

%     speed = abs((diff(colPlot)+ diff(rowPlot) + diff(zPlot))/1000/0.005);
    speedz = abs(diff(zPlot))/1000/0.005;
    %speedR = diff(rPlot)/1000/0.005;


    allSpeed = [allSpeed; speedz(:)];
    
end