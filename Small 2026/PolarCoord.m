
%Spherical coordinate
%% vector plot
threshold =20;
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
%% binning for 2D vector plot XY
%==> this would be better to do in spherical coordinate
binX = linspace(round(min(filterPos(:,4))),round(max(filterPos(:,4))),40);
binZ = linspace(round(min(filterPos(:,5))),round(max(filterPos(:,5))),40);

[X,Z] = meshgrid(binX(1:end-1),binZ(1:end-1));

U =zeros(size(X));
V = zeros(size(Z));
for i =1:length(binX)-1
    for j =1:length(binZ)-1

        cIdx = and(filterPos(:,4)>binX(i), filterPos(:,4)<=binX(i+1));
        cIdz  = and(filterPos(:,5)>binZ(j), filterPos(:,5)<=binZ(j+1));
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

% figure