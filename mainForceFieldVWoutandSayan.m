%% Force fields plotting

%% input
InputDirectory = 'F:\WoutVT\20251003_trapping_statisics_BE_only_2 (1 MHz 2V)\200nm_PS_1MHz_2V';

load([InputDirectory '\' 'trackResults.mat']);

calc = 'speed'; %'speed','accForce' or 'dragForce' 
TrackedData = trackRes.traces;
%color = getColorFromCmap('lines',length(TrackedData));
expTime = 0.01; %in s
density = 1055;% in kg/m3
R = 125*10^-9; %hydrodynamic in meter
volume = 4/3*pi*R^3;%in m3
mass = volume*density;%in kg
viscosity = 0.001;% in Pa.s

%% section to filter a specific subset of particles (e.g. high directionality)
% %% Filter out traces that are not going to the trap
for i=1:length(TrackedData)
    
    currTrace = TrackedData{i,1};
    % test1 = ~all(currTrace.t<40);
    test2 = height(currTrace)>5;

    %test= sum(gradient(currTrace.z));
%     test1 = sum(gradient(currTrace.z))>2000;
%     test2 = norm([currTrace.col(1), currTrace.row(1)]-CM)> norm([currTrace.col(end), currTrace.row(end)]-CM);
%     test3 = height(currTrace)>5;
%     test4 = and(currTrace.col(end)< CM(1)+2000,currTrace.col(end)> CM(1)-2000); 
%     test5 = and(currTrace.row(end)< CM(2)+2000,currTrace.row(end)> CM(2)-2000); 
%     test = logical(test1*test2*test3*test4*test5);
    test = logical(test2);
    
    if test
%         figure(1)
%         hold on
%         
%         plot3(currTrace.col,currTrace.row,currTrace.z);
%         axis image
    else
        TrackedData{i,1} = [];
    end
    
    
end
%%
traces2Keep = TrackedData(~cellfun(@isempty,TrackedData(:,1)),1);
FullPos = [];
figure
hold on
for i=1:size(traces2Keep,1)
    track = traces2Keep{i,1};
    
    switch calc
        case 'speed'
            %calculate speed from first derivative of displacements
            ax = diff(track.col)/1000/expTime;
            ay = diff(track.row)/1000/expTime;
            az = diff(track.z)/1000/expTime;
        case 'accForce'    
            %calculate force based on acceleration ==> Second derivative of
            %displacement
            ax = diff(diff(track.col))*10^-9/(expTime^2)*mass;
            ay = diff(diff(track.row))*10^-9/(expTime^2)*mass;
            az = diff(diff(track.z))*10^-9/(expTime^2)*mass;
        case 'dragForce'
        %calculate force based drag force
            ax = 6 * pi * viscosity * R* diff(track.col)*10^-9/(expTime);
            ay = 6 * pi * viscosity * R* diff(track.row)*10^-9/(expTime);
            az = 6 * pi * viscosity * R* diff(track.z)*10^-9/(expTime);
        otherwise
            error('Unexpected metric requested for calculation, only accept speed accForce and dragForce')
    end
    
    x = track.col;
    y = track.row;
    z = track.z;
    %plot(track.z(2:end),ax+ay+az)
    FullPos = [FullPos; [ax,ay,az , x(2:end), y(2:end),z(2:end)]];
    
end



%% get max Force

normForce = sqrt(FullPos(:,1).^2+FullPos(:,2).^2+FullPos(:,3).^2)./10^-12;
normSpeed = sqrt(FullPos(:,1).^2+FullPos(:,2).^2+FullPos(:,3).^2);




%% tmp tr
figure
quiver3(FullPos(:,4),FullPos(:,5),FullPos(:,6),FullPos(:,1),FullPos(:,2),FullPos(:,3))
axis image
box on
set(gcf,'color','w')

%%

close all

cst = 10;
figure
quiverC3D(FullPos(:,4),FullPos(:,5),FullPos(:,6),FullPos(:,1)*cst,FullPos(:,2)*cst,FullPos(:,3)*cst)
axis image
view(3)
colorbar
colormap('jet')
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
