%% Force fields plotting

%% input
InputDirectory = 'F:\Colloidal memory\2024-10-25\Sample 3\Trigger';

load([InputDirectory '\' 'trackResults.mat']);

calc = 'speed'; %'speed','accForce' or 'dragForce' 
TrackedData = trackRes.traces;
%color = getColorFromCmap('lines',length(TrackedData));
expTime = 0.01; %in s
density = 1055;% in kg/m3
R = 40*10^-9; %hydrodynamic in meter
volume = 4/3*pi*R^3;%in m3
mass = volume*density;%in kg
viscosity = 0.001;% in Pa.s

%% section to filter a specific subset of particles (e.g. high directionality)
% %% Filter out traces that are not going to the trap
for i=1:length(TrackedData)
    
    currTrace = TrackedData{i,1};
    test1 = ~all(currTrace.t<40);
    test2 = height(currTrace)>5;

    %test= sum(gradient(currTrace.z));
%     test1 = sum(gradient(currTrace.z))>2000;
%     test2 = norm([currTrace.col(1), currTrace.row(1)]-CM)> norm([currTrace.col(end), currTrace.row(end)]-CM);
%     test3 = height(currTrace)>5;
%     test4 = and(currTrace.col(end)< CM(1)+2000,currTrace.col(end)> CM(1)-2000); 
%     test5 = and(currTrace.row(end)< CM(2)+2000,currTrace.row(end)> CM(2)-2000); 
%     test = logical(test1*test2*test3*test4*test5);
    test = logical(test1*test2);
    
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
            ax = diff(diff(track.row))*10^-9/(expTime^2)*mass;
            ay = diff(diff(track.col))*10^-9/(expTime^2)*mass;
            az = diff(diff(track.z))*10^-9/(expTime^2)*mass;
        case 'dragForce'
        %calculate force based drag force
            ax = 6 * pi * viscosity * R* diff(track.col)*10^-9/(expTime);
            ay = 6 * pi * viscosity * R* diff(track.row)*10^-9/(expTime);
            az = 6 * pi * viscosity * R* diff(track.z)*10^-9/(expTime);
        otherwise
            error('Unexpected metric requested for calculation, only accept speed accForce and dragForce')
    end
    
    x = track.row;
    y = track.col;
    z = track.z;
    plot(track.z(2:end),ax+ay+az)
    FullPos = [FullPos; [ax,ay,az , x(2:end), y(2:end),z(2:end)]];
    
end
%% tmp tr
max(FullPos(:,6));
quiver3(FullPos(:,4),FullPos(:,5),FullPos(:,6),FullPos(:,1),FullPos(:,2),FullPos(:,3))
axis image

%% 
cst = 10;
close all
colormap jet
quiverC3D(FullPos(:,4),FullPos(:,5),FullPos(:,6),FullPos(:,1)*cst,FullPos(:,2)*cst,FullPos(:,3)*cst)
axis equal
%%

step = 200;
full = [];
for x = min(FullPos(:,4)):step:max(FullPos(:,4))
    for y = min(FullPos(:,5)):step:max(FullPos(:,5))
        for z = 0:-step:-4000
            indx = find(FullPos(:,4)>x & FullPos(:,4)<x+step);
            indy = find(FullPos(:,5)>y & FullPos(:,5)<y+step);
            indz = find(FullPos(:,6)>z & FullPos(:,6)<z+step);
            indxy = intersect(indx,indy);
            indxyz = intersect(indxy,indz);
             if ~isempty(indxyz) 
                 full = [full ;x,y,z,mean(FullPos(indxyz,1)),mean(FullPos(indxyz,2)),mean(FullPos(indxyz,3))];
             end
        end
    end
end

% %clean data
% ind = find(or(full(:,6)<-3e-10,full(:,6)>5e-10));
% full(ind,:)= [];
% ind = find(or(full(:,5)<-3e-10,full(:,5)>5e-10));
% full(ind,:) = [];
% ind = find(or(full(:,4)<-3e-10,full(:,4)>5e-10));
% full(ind,:) = [];

close all
colormap jet
quiverC3D(full(:,1),full(:,2),full(:,3),full(:,4)*2*10^15,full(:,5)*2*10^15,full(:,6)*2*10^15)
axis equal
% xlim([6000 13000])
% ylim([20000 27000])




%%
test = unique(full(:,3));
avgZ = zeros(size(test));
for i = 1:length(test)
avgZ(i) = nanmean(full(full(:,3)==test(i),6));
end
figure
plot(test,avgZ);
xlabel('z position')
ylabel('Average Force')
