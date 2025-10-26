%% Force fields plotting

%% input
% InputDirectory = 'F:\Colloidal memory\2024-10-25\Sample 3\Trigger';
% 
% load([InputDirectory '\' 'trackResults.mat']);

calc = 'speed'; %'speed','accForce' or 'dragForce' 
TrackedData = trackRes.traces;
%color = getColorFromCmap('lines',length(TrackedData));
expTime = 0.008; %in s
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
% figure
% hold on
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


%% Rotate data in XY
angle =45;%45 for DEP in quadrupole
clear rotated;
threshold = 50;
%RotMat = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];

RotMat = rotz(angle);

rotated(:,1:3) =FullPos(:,1:3) *RotMat;

rotated(:,4:6) = FullPos(:,4:6) *RotMat;

cst = 10;
figure
quiverC3D(rotated(:,4),rotated(:,5), rotated(:,6),rotated(:,1)*cst,rotated(:,2)*cst,rotated(:,3)*cst)
axis image
view(2)
colorbar
colormap('jet')
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image


%24000 32000 pDEP % 26500 31000 nDEP
rotated(or(rotated(:,4)<24500,rotated(:,4)>33000),:) = [];

%sum of XYZ should be >100
rotated(sum(abs(rotated(:,1:3)),2)<threshold,:)=[];


cst = 2;
figure
quiverC3D(rotated(:,4),rotated(:,6),zeros(size(rotated(:,4))),rotated(:,1)*cst,rotated(:,3)*cst,zeros(size(rotated(:,4)))*cst,200)
axis image
view(2)
colorbar
colormap('jet')
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image

%% histogram of rotated

binX = linspace(round(min(rotated(:,4))),round(max(rotated(:,4))),40);
binZ = linspace(round(min(rotated(:,6))),round(max(rotated(:,6))),40);

[X,Z] = meshgrid(binX(1:end-1),binZ(1:end-1));

U =zeros(size(X));
V = zeros(size(Z));
for i =1:length(binX)-1
    for j =1:length(binZ)-1
        cIdx = and(rotated(:,4)>binX(i), rotated(:,4)<=binX(i+1));
        cIdz  = and(rotated(:,6)>binZ(j), rotated(:,6)<=binZ(j+1));
        cId = and(cIdx,cIdz);
        U(j,i) = mean(rotated(cId,1));
        V(j,i) = mean(rotated(cId,3));
        
    end
end
% 
% figure
% quiver(X,Z,U,V)
cst = 5;

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
xlim([24000,32000])

%%
figure
streamline(X,Z,U,V,X,Z)

%% 
  
hold on

quiverC3D(X,Z,zeros(size(X)),U*cst,V*cst,zeros(size(Z))*cst)
axis image
view(2)
colorbar
colormap('jet')

xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
xlim([24000,32000])

%% rotation along the electrode

angle =45;%45 for DEP 
clear rotated;
threshold = 50;
%RotMat = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];

RotMat = rotz(angle);

rotated(:,1:3) =FullPos(:,1:3) *RotMat;

rotated(:,4:6) = FullPos(:,4:6) *RotMat;


%24000 32000 pDEP % 26500 31000 nDEP
rotated(or(rotated(:,5)<4000,rotated(:,5)>8500),:) = [];

%sum of XYZ should be >100
rotated(sum(abs(rotated(:,1:3)),2)<threshold,:)=[];

binX = linspace(round(min(rotated(:,4))),round(max(rotated(:,4))),40);
binZ = linspace(round(min(rotated(:,6))),round(max(rotated(:,6))),40);

[X,Z] = meshgrid(binX(1:end-1),binZ(1:end-1));

U =zeros(size(X));
V = zeros(size(Z));
for i =1:length(binX)-1
    for j =1:length(binZ)-1

        cIdx = and(rotated(:,4)>binX(i), rotated(:,4)<=binX(i+1));
        cIdz  = and(rotated(:,6)>binZ(j), rotated(:,6)<=binZ(j+1));
        cId = and(cIdx,cIdz);
        U(j,i) = mean(rotated(cId,1));
        V(j,i) = mean(rotated(cId,3));
        
    end
end
% 
% figure
% quiver(X,Z,U,V)
cst = 3;

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








%% Top view

CM = [24400 16400];
threshold = 1e-13;
% 
ind2Plot =[abs(FullPos(:,1))> threshold, abs(FullPos(:,2))> threshold, abs(FullPos(:,3))> threshold];

ind2Plot = (sum(ind2Plot,2))>0;

figure
cst = 25;
close all
colormap jet
quiverC3D(FullPos(ind2Plot,4),FullPos(ind2Plot,5),FullPos(ind2Plot,6),FullPos(ind2Plot,1)*cst,FullPos(ind2Plot,2)*cst,FullPos(ind2Plot,3)*cst)
axis image
ylim([8400 39400])
xlim([800 32000])
view(2)
box on
colorbar
set(gcf,'color','w')
%% 3D View
figure
colormap jet
quiverC3D(FullPos(ind2Plot,4),FullPos(ind2Plot,5),FullPos(ind2Plot,6),FullPos(ind2Plot,1)*cst,FullPos(ind2Plot,2)*cst,FullPos(ind2Plot,3)*cst)
axis image

view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
r =5000;
xlim([CM(2)-r,CM(2)+r])
ylim([CM(1)-r,CM(1)+r])
view(3)
box on
set(gcf,'color','w')


%% 
figure
colormap jet
quiverC3D(FullPos(ind2Plot,4),FullPos(ind2Plot,5),FullPos(ind2Plot,6),FullPos(ind2Plot,1)*cst,FullPos(ind2Plot,2)*cst,FullPos(ind2Plot,3)*cst)
axis image

view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image
r = 3000;
xlim([CM(2)-10000-r,CM(2)-10000+r])
ylim([CM(1)-10000-r,CM(1)-10000+r])
zlim([-2200 2300])
view(135.8,-0)
%view(gca,[90 0])
box on
set(gcf,'color','w')
















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
