%% Load data
clear; clc;

%%%%user input
InputDirectory = 'C:\Users\sayan\OneDrive\Desktop\KU Leuven\data\200nm\data1\100kHz 5V';
load([InputDirectory '\' 'trackResults.mat']);
    beep on; beep; 
%% User inputs
minpoints=10;         %minimum no. of points to be there in each trajectory
zstart = 1700;        %particle is atleast below this height at the beginning  
%zcond=1200;          %minimum z to consider traped, 20% points are above this height  %%decrease to get more traces
t_detrap=20;          %time when detrap starts
calc = 'speed';       %'speed','accForce' or 'dragForce' 
expTime = 10;         %in ms
R = 125*10^-9;        %hydrodynamic radius in meter
density = 1055;       % in kg/m3

volume = 4/3*pi*R^3;%in m3
mass = volume*density;%in kg
viscosity = 0.001;% in Pa.s
expTime=expTime/1000;% in s
num_videos=trackRes.traces{end,2};
%% Filter out the good trajectories

TrackedData = trackRes.traces;
for i=1:length(TrackedData)

    currTrace = TrackedData{i,1};          % go through each of the x12 arrays. 
    test1=height(currTrace)>minpoints;     % the trace has atleast minpoints data points
    %test2=1;
    test2=currTrace.rT(end)/1000<=t_detrap; % the end of the trace is before detrapping starts
    test3=1;
    %test3=currTrace.z(1)<=zstart;          % particle is not already trapped
    test4=1;
    %test4=currTrace.z(end)>=zcond;         % final z is above zcond, particle remains trapped
    test5=1;
    %test5=(length(find(currTrace.z>zcond)))>=0.50*height(currTrace); %atleast 20 percent of the trace is above zcond
    
    test=logical(test1*test2*test3*test4*test5);

    if ~test
        TrackedData{i,1} = [];
    end

    if(currTrace.rT(1)/1000<t_detrap && currTrace.rT(end)/1000>t_detrap)
        idx=find(currTrace.rT/1000 >= 20, 1, 'first');
        currTrace=currTrace(1:idx,:);
        if(height(currTrace)>minpoints)
        TrackedData{i,1}=currTrace;
        end
    end

    if(currTrace.rT(1)/1000>t_detrap)
        TrackedData{i,1} = [];
    end

 end

GoodPos = TrackedData(~cellfun(@isempty,TrackedData(:,1)),1);

%% Colour z conc Flip idea

nData =sum(cellfun(@height,GoodPos));
AllPos=zeros(nData,4);
for i=1:length(GoodPos)
    currTrace=GoodPos{i,1};
    p=find(AllPos==0,1);
    AllPos(p:p+height(currTrace)-1,1)=currTrace.row;
    AllPos(p:p+height(currTrace)-1,2)=currTrace.col;
    AllPos(p:p+height(currTrace)-1,3)=currTrace.z;
    AllPos(p:p+height(currTrace)-1,4)=currTrace.rT;
end

points_per_sec= 50;
AllPos(:,4) = round(AllPos(:,4)./1000, ceil(log10(points_per_sec)));

if(floor(min(AllPos(:,3)))<-2000)
    z_start=floor(min(AllPos(:,3)));
else
    z_start =-2000;
end

if(ceil(max(AllPos(:,3)))>2000)
    z_end=ceil(max(AllPos(:,3)));
else
    z_end =2000;
end



z_div_num=z_end-z_start+1;
max_time = max(AllPos(:,4));
time_div_num = round(max_time*points_per_sec)+1;

z_t=zeros(z_div_num,time_div_num);

tic
for i = 1 : length(AllPos)
    m=z_end-round(AllPos(i,3))+1;
    n=round((AllPos(i,4))*points_per_sec)  +1;

    z_t(m,n)=z_t(m,n)+1;

end
toc



time=0:(1/points_per_sec):max_time;
%time=(time(:)*30/max_time)';
z=z_start:1:z_end;
z_t=flipud(z_t);
 

figure(1)
imagesc(time, z/1000,z_t)
clim([0 5])
colormap turbo
colorbar
axis xy;
xlabel('Time (s)','FontWeight','bold','FontSize',16)
ylabel('Z (\mum)','FontWeight','bold','FontSize',16)
%title('60nm 15MHz ','FontWeight','bold','FontSize',16)
set(gca,'fontweight','bold','LineWidth',2.5,'Box','on','XColor',[0,0,0],'YColor',[0,0,0]);
ax = gca;
ax.FontSize = 14; 
ax.FontName='Times New Roman';
%% Z Conc Pablo idea

nData =sum(cellfun(@height,GoodPos));
AllPos=zeros(nData,4);
for i=1:length(GoodPos)
    currTrace=GoodPos{i,1};
    p=find(AllPos==0,1);
    AllPos(p:p+height(currTrace)-1,1)=currTrace.row;
    AllPos(p:p+height(currTrace)-1,2)=currTrace.col;
    AllPos(p:p+height(currTrace)-1,3)=currTrace.z;
    AllPos(p:p+height(currTrace)-1,4)=currTrace.rT;
end


num_lines=8;
min_z=-2000;
max_z=2000;
sampling=5; % number of time points per second

z_breaks=linspace(min_z,max_z,num_lines+1);

time=0:1/sampling:round(max(AllPos(:,4))/1000,ceil(log10(sampling)));

z_conc=cell(num_lines,1);

for i=1:height(z_conc)
    z_conc{i}=zeros(1,length(time));
end

for i=1:length(AllPos)
  
    for j=1:length(z_breaks)-1

    if(AllPos(i,3)>z_breaks(j)&&AllPos(i,3)<=z_breaks(j+1))
        
        t=round((AllPos(i,4)/1000),ceil(log10(sampling)));
        z_conc{j}(1,fix((sampling*t)+1))=z_conc{j}(1,fix((sampling*t)+1))+1;

    end

    end
    
end

for i=1:height(z_conc)
    z_conc{i}(1:length(time))=z_conc{i}(1:length(time))./num_videos;
end


legendText=cell(1,num_lines);
colors = brighten(jet(num_lines),0); 
figure
hold on
for i=1:num_lines
semilogy(time,z_conc{i},'-','LineWidth',2,'Color', colors(i,:))
legendText{i}=sprintf(' %.2f to %.2f', z_breaks(i),z_breaks(i+1));
end
hold off
xlabel('Time (s)','FontSize',14,'FontWeight','bold')
ylabel('Number','FontSize',14,'FontWeight','bold')
legend(legendText, 'Location', 'northeastoutside',FontWeight='bold',FontSize=9,EdgeColor='none');
set(gca,'fontweight','bold','FontSize',14,'LineWidth',2,'Box','on','XColor',[0,0,0],'YColor',[0,0,0]);




for i=1:num_lines
    figure
plot(time,z_conc{i},'-','LineWidth',2,'Color', colors(i,:))
title(sprintf(' %.2f to %.2f', z_breaks(i),z_breaks(i+1)),'FontWeight','bold')
xlabel('Time (s)')
ylabel('Number')
end


%% Mask create

% nData =sum(cellfun(@height,GoodPos));
% AllPos=zeros(nData,4);
% for i=1:length(GoodPos)
%     currTrace=GoodPos{i,1};
%     p=find(AllPos==0,1);
%     AllPos(p:p+height(currTrace)-1,1)=currTrace.row;
%     AllPos(p:p+height(currTrace)-1,2)=currTrace.col;
%     AllPos(p:p+height(currTrace)-1,3)=currTrace.z;
%     AllPos(p:p+height(currTrace)-1,4)=currTrace.rT;
% end
% 
% figure
% histogram(AllPos(:,1),100);
% figure
% histogram(AllPos(:,2),100);
% figure
% histogram(AllPos(:,3),100);
% 
% figure
% [N, edges] = histcounts(AllPos(:,1), 100);
% x = edges(1:end-1) + diff(edges)/2; 
% y = N;
% dmin = 3500; % separation between Gaussians (set to slightly less than expected spacing)
% 
% [pks, locsx] = findpeaks(y, x, 'MinPeakDistance', dmin);
% plot(x, y, 'k-', 'LineWidth', 2)
% hold on
% plot(locsx, pks, 'ro', 'MarkerSize', 10, 'LineWidth', 2)
% xlabel('X Value')
% ylabel('Counts')
% hold on
% window = 1000;   % half-width around each peak (adjust!)
% 
% for i = 1:length(locsx)
% 
%     idx = abs(x - locsx(i)) < window;
% 
%     % Gaussian fit
%     f = fit(x(idx)', y(idx)', 'gauss1');
% 
%     % Smooth x for plotting the fit
%     xf = linspace(min(x(idx)), max(x(idx)), 50000);
%     yf = f(xf);
% 
%     % Plot Gaussian
%     plot(xf, yf,'b', 'LineWidth', 1)
% 
%     % (optional) store parameters
%     mux(i)    = f.b1;
%     sigmax(i) = f.c1/sqrt(2);
%     Ax(i)     = f.a1;
% 
%     % Gaussian center
%     mu = f.b1;
% 
%     % Value of Gaussian at center (for correct y-position)
%     ymu = f(mu);
% 
%     % Plot center as blue circle
%     plot(mu, ymu, 'bo', 'MarkerSize', 8, 'LineWidth', 2)
% end
% hold on
% xline(mux-2*sigmax);
% hold on
% xline(mux+2*sigmax);
% 
% 
% 
% figure
% [N, edges] = histcounts(AllPos(:,2), 100);
% x = edges(1:end-1) + diff(edges)/2; 
% y = N;
% dmin = 3500; % separation between Gaussians (set to slightly less than expected spacing)
% 
% [pks, locsy] = findpeaks(y, x, 'MinPeakDistance', dmin);
% plot(x, y, 'k-', 'LineWidth', 2)
% hold on
% plot(locsy, pks, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
% xlabel('Y Value')
% ylabel('Counts')
% hold on
% window = 500;   % half-width around each peak (adjust!)
% 
% for i = 1:length(locsy)
% 
%     idx = abs(x - locsy(i)) < window;
% 
%     % Gaussian fit
%     f = fit(x(idx)', y(idx)', 'gauss1');
% 
%     % Smooth x for plotting the fit
%     xf = linspace(min(x(idx)), max(x(idx)), 500000);
%     yf = f(xf);
% 
%     % Plot Gaussian
%     plot(xf, yf,'b', 'LineWidth', 1)
% 
%     % (optional) store parameters
%     muy(i)    = f.b1;
%     sigmay(i) = f.c1/sqrt(2);
%     Ay(i)     = f.a1;
% 
%     % Gaussian center
%     mu = f.b1;
% 
%     % Value of Gaussian at center (for correct y-position)
%     ymu = f(mu);
% 
%     % Plot center as blue circle
%     plot(mu, ymu, 'bo', 'MarkerSize', 8, 'LineWidth', 2)
% end
% hold on
% xline(muy-2*sigmay);
% hold on
% xline(muy+2*sigmay);
% 
% 
% 
% figure
% [N, edges] = histcounts(AllPos(:,3), 100);
% x = edges(1:end-1) + diff(edges)/2; 
% y = N;
% dmin = 500; % separation between Gaussians (set to slightly less than expected spacing)
% 
% [pks, locsz] = findpeaks(y, x, 'MinPeakDistance', dmin,'MinPeakHeight',length(AllPos)/500);
% plot(x, y, 'k-', 'LineWidth', 2)
% hold on
% plot(locsz, pks, 'ro', 'MarkerSize', 10, 'LineWidth', 2)
% xlabel('X Value')
% ylabel('Counts')
% hold on
% window = 200;   % half-width around each peak (adjust!)
% 
% for i = 1:length(locsz)
% 
%     idx = abs(x - locsz(i)) < window;
% 
%     % Gaussian fit
%     f = fit(x(idx)', y(idx)', 'gauss1');
% 
%     % Smooth x for plotting the fit
%     xf = linspace(min(x(idx)), max(x(idx)), 50000);
%     yf = f(xf);
% 
%     % Plot Gaussian
%     plot(xf, yf,'b', 'LineWidth', 1)
% 
%     % (optional) store parameters
%     muz(i)    = f.b1;
%     sigmaz(i) = f.c1/sqrt(2);
%     Az(i)     = f.a1;
% 
%     % Gaussian center
%     mu = f.b1;
% 
%     % Value of Gaussian at center (for correct y-position)
%     ymu = f(mu);
% 
%     % Plot center as blue circle
%     plot(mu, ymu, 'bo', 'MarkerSize', 8, 'LineWidth', 2)
% end
% hold on
% xline(muz-3*sigmaz);
% hold on
% xline(muz+3*sigmaz);
% 
% Mask1 = zeros(round( max(AllPos(:,1))/1 ), round( max(AllPos(:,2)/1 )));
% 
% for i=1:length(locsx)
% 
%     for j=1:length(locsy)
%         xmin=round( (mux(i)-1.5*sigmax(i))/1 );
%         xmax=round( (mux(i)+1.5*sigmax(i))/1 );
%         ymin=round( (muy(j)-1.5*sigmay(j))/1 );
%         ymax=round( (muy(j)+1.5*sigmay(j))/1 );
%         Mask1(xmin:xmax,ymin:ymax)=1;
%     end
% 
% end
% figure
% imagesc(Mask1( 1:round(length(Mask1)/2) , 1:round(width(Mask1)/2) ));
% axis image
% height_filter=muz(end)-2*sigmaz(end);

%% calculate speed/ force (not used now, but need to be run for next section to work)

VelPostime = cell(length(GoodPos),1);
accFPostime = cell(length(GoodPos),1);
dragFPostime = cell(length(GoodPos)-1,1);

Vel = [];
accForce=[];
dragForce=[];

for i=1:size(GoodPos,1)
    currTrace = GoodPos{i,1};
    currX=currTrace.col;
    currY=currTrace.row;
    currZ=currTrace.z;
    currT=currTrace.rT/1000;    %rT is in ms, convert it to s by rT/1000

    switch calc
        case 'speed'
            vx=diff(currX)./diff(currT);
            vy=diff(currY)./diff(currT);
            vz=diff(currZ)./diff(currT);
            
            Vel = [Vel; [vx,vy,vz , currX(1:end-1)+diff(currX)/2, currY(1:end-1)+diff(currY)/2,currZ(1:end-1)+diff(currZ)/2, currT(1:end-1)]];
            varNames = ["vx","vy","vz","v_mag","x","y","z","time"];
            currVel = table(vx,vy,vz,sqrt(vx.^2 + vy.^2 + vz.^2) ,currX(1:end-1)+diff(currX)/2, currY(1:end-1)+diff(currY)/2, currZ(1:end-1)+diff(currZ)/2 ,currT(1:end-1),'VariableNames',varNames);
            VelPostime{i,1}=currVel;

        case 'accForce'
            vx=diff(currX)./diff(currT); 
            vy=diff(currY)./diff(currT); 
            vz=diff(currZ)./diff(currT); 
            t_vel=diff(currT);

            fx=mass*(diff(vx)./diff(t_vel));
            fy=mass*(diff(vx)./diff(t_vel));
            fz=mass*(diff(vx)./diff(t_vel));

            xpos_vel=currX(1:end-1)+diff(currX)/2;
            ypos_vel=currY(1:end-1)+diff(currY)/2;
            zpos_vel=currZ(1:end-1)+diff(currZ)/2;

            accForce=[accForce; [fx,fy,fz, xpos_vel(1:end-1)+diff(xpos_vel)/2, ypos_vel(1:end-1)+diff(ypos_vel)/2, zpos_vel(1:end-1)+diff(zpos_vel)/2, currT(1:end-2)]];
            
            varNames = ["fx","fy","fz","f_mag","x","y","z","time"];
            currF = table(fx,fy,fz, sqrt(fx.^2 + fy.^2 + fz.^2) , xpos_vel(1:end-1)+diff(xpos_vel)/2, ypos_vel(1:end-1)+diff(ypos_vel)/2, zpos_vel(1:end-1)+diff(zpos_vel)/2 ,currT(1:end-2),'VariableNames',varNames);
            accFPostime{i,1}=currF;

        case 'dragFroce'
            vx=diff(currX)./diff(currT);
            vy=diff(currY)./diff(currT);
            vz=diff(currZ)./diff(currT);

            fx=6*pi*viscosity*R*vx;
            fy=6*pi*viscosity*R*vy;
            fz=6*pi*viscosity*R*vz;

            dragForce=[dragFroce; [fx,fy,fz, currX(1:end-1)+diff(currX)/2, currY(1:end-1)+diff(currY)/2,currZ(1:end-1)+diff(currZ)/2, currT(1:end-1)]]; 
            varNames = ["fx","fy","fz","f_mag" ,"x","y","z","time"];
            currF = table(fx,fy,fz,sqrt(fx.^2 + fy.^2 + fz.^2), currX(1:end-1)+diff(currX)/2, currY(1:end-1)+diff(currY)/2, currZ(1:end-1)+diff(currZ)/2 ,currT(1:end-1),'VariableNames',varNames);
            dragFPostime{i,1}=currF;
        otherwise
            error('Wrong option, choose from speed accForce and dragForce');
    end
end
%% Position plotting

std_x=zeros(1,length(GoodPos));
std_y=zeros(1,length(GoodPos));
std_z=zeros(1,length(GoodPos));
mean_z=zeros(1,length(GoodPos));

x_position=zeros(1,length(GoodPos));
y_position=zeros(1,length(GoodPos));

max_x=zeros(1,length(GoodPos));
min_x=zeros(1,length(GoodPos));
max_y=zeros(1,length(GoodPos));
min_y=zeros(1,length(GoodPos));


div=1; % no. of divisions to the data u want to make (leave at 1 for whole data)
s=1; %1 to div
range=[(round(length(GoodPos)/div)*(s-1))+1 round(length(GoodPos)/div)*s];
trajectories=cell(range(2),1);

for i=range(1):range(2)
    currTrace = GoodPos{i,1};
    currVel = VelPostime{i,1};
    currTrace.rT=currTrace.rT/1000;
    x1=currTrace.row - currTrace.row(end);
    y1=currTrace.col - currTrace.col(end);
    z1=currTrace.z - currTrace.z(end);
    t1=currTrace.rT - currTrace.rT(1)+0.01;
    t1=t1./t1(end);

   
    
    x=currTrace.row;
    y=currTrace.col;
    z=currTrace.z;
    r_vec=sqrt(x1.^2+y1.^2+z1.^2);
    t=currTrace.rT;

    max_x(i)=max(x);
    max_y(i)=max(y);

    min_x(i)=min(x);
    min_y(i)=min(y);

    std_x(i)=std(x1);
    std_y(i)=std(y1);
    std_z(i)=std(z);
    mean_z(i)=mean(z);

    trajectories{i} = table(t(:), x1(:), y1(:), z(:), 'VariableNames', {'t','x','y','z'});

    x_position(i)=mean( x(end-round(length(x)/2):end) );
    y_position(i)=mean( y(end-round(length(y)/2):end) );

    %currVel.time=currVel.time-currVel.time(1);
    %currVel.v_mag=currVel.v_mag-currVel.v_mag(end);
    vx=(currVel.vx);
    vy=(currVel.vy);
    vz=(currVel.vz);
    vmag=currVel.v_mag;
    v_t=currVel.time;
    v_t1=currVel.time-currVel.time(1);
    v_t1=v_t1./v_t1(end);

    vx_pos=currVel.x-currVel.x(end);
    vy_pos=currVel.y-currVel.y(end);
    vz_pos=currVel.z;

    [theta, rho, z_cyl] = cart2pol(x1, y1, z);
    % rho = sqrt(x1.^2 + y1.^2 + z.^2);
    % theta = atan2(y1, x1);
    % phi = acos(z./rho);
    % [x_cart, y_cart, z_cart] = pol2cart(theta, rho, z_cyl);


   %%%% Plotting start

    % figure(20)
    % polarplot(theta, rho);
    % patch([theta' nan],[rho' nan], [t1' nan],'EdgeColor','interp','FaceColor','none')
    % hold on

%     figure(1)
%     plot3(x,y,z);
%     patch([x' nan],[y' nan],[z' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
%     view(3)
%     xlabel('X (nm)','FontSize', 15, 'FontWeight', 'bold');
%     ylabel('Y (nm)','FontSize', 15, 'FontWeight', 'bold');
%     zlabel('Z (nm)','FontSize', 15, 'FontWeight', 'bold');
%     set(gca,'fontweight','bold','LineWidth',2.5,'Box','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0]);
%     ax = gca;
%     ax.FontSize = 20; 
%     ax.FontName='Times New Roman';
%     colorbar("east");
%     title('3D plot of position (not centred)')
%     hold on

%     figure(15)
%     plot3(x1,y1,z);
%     patch([x1' nan],[y1' nan],[z' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
%     view(3)
%     xlabel('X (nm)','FontSize', 15, 'FontWeight', 'bold');
%     ylabel('Y (nm)','FontSize', 15, 'FontWeight', 'bold');
%     zlabel('Z (nm)','FontSize', 15, 'FontWeight', 'bold');
%     set(gca,'fontweight','bold','LineWidth',2.5,'Box','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0]);
%     ax = gca;
%     ax.FontSize = 20; 
%     ax.FontName='Times New Roman';
%     colorbar("east");
%     title('3D plot of position (centred)')
%     hold on


%    % 
%     figure(2)
%     plot(x1,z) 
%     patch([x1' nan],[z' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
%     xlabel('X (\mum)','FontSize', 15, 'FontWeight', 'bold');
%     ylabel('Z (\mum)','FontSize', 15, 'FontWeight', 'bold');
%     xticks([ -10000 -5000 0 5000 10000 15000 20000])
%     xticklabels({'-10','-5','0','5','10', '15','20'})
%     yticks([-2000 -1500 -1000 -500 0 500 1000 1500 2000])
%     yticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'})
%     set(gca,'fontweight','bold','LineWidth',2,'Box','on','XColor',[0,0,0],'YColor',[0,0,0]);
%     ax = gca;
%     ax.FontSize = 14; 
%     ax.FontName='Times New Roman';
%     title('X-Z Plot')
%     c = colorbar(ax);
% %ax.Visible = 'off';
% c.FontWeight="bold";
% c.FontName='Times New Roman';
% c.FontSize=12;
% c.LineWidth = 1;
% c.Label.String = 'Relative  Time';
% c.Label.FontSize = 16;
%     hold on
% 
%     figure(3)
%     plot(y1,z) %t1 = all at beginning random. t = random lines join as they get trapped
%     patch([y1' nan],[z' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
%     xlabel('Y (\mum)','FontSize', 15, 'FontWeight', 'bold');
%     ylabel('Z (\mum)','FontSize', 15, 'FontWeight', 'bold');
%     xticks([ -10000 -5000 0 5000 10000 15000 20000])
%     xticklabels({'-10','-5','0','5','10', '15','20'})
%     yticks([-2000 -1500 -1000 -500 0 500 1000 1500 2000])
%     yticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'})
%     set(gca,'fontweight','bold','LineWidth',2,'Box','on','XColor',[0,0,0],'YColor',[0,0,0]);
%     ax = gca;
%     ax.FontSize = 14; 
%     ax.FontName='Times New Roman';
%     title('Y-Z Plot')
%     c = colorbar(ax);
% %ax.Visible = 'off';
% c.FontWeight="bold";
% c.FontName='Times New Roman';
% c.FontSize=12;
% c.LineWidth = 1;
% c.Label.String = 'Relative  Time';
% c.Label.FontSize = 16;
%     hold on
% 
%    % 
%     figure(4)
%     plot(x,y) %t =blue also in center, t1= all blue outside, yellow inside
%     patch([x' nan],[y' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
%     xlabel('X (nm)','FontSize', 15, 'FontWeight', 'bold');
%     ylabel('Y (nm)','FontSize', 15, 'FontWeight', 'bold');
%     set(gca,'fontweight','bold','LineWidth',2.5,'Box','off','XColor',[0,0,0],'YColor',[0,0,0]);
%     ax = gca;
%     ax.FontSize = 20; 
%     ax.FontName='Times New Roman';
%     title('X-Y Top view Plot')
%     hold on
% 
%     figure(5)
%     plot(t,z) %t=random lines join as they get trapped, t1= all random at beginning
%     patch([t' nan],[z' nan],[t' nan],'EdgeColor','interp','FaceColor','none')
%     title('Z vs Time plot');
%     hold on




    % figure(6)
    % plot(t1,z) %t=random lines join as they get trapped, t1= all random at beginning
    % patch([t1' nan],[z' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
    % title('Z vs Time plot');
    % hold on
   % 
   %  figure(6)
   %  plot(t1,r_vec);%t=random lines join as they get trapped, t1= all random at beginning
   %  patch([t1' nan],[r_vec' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
   %  title('r vs Time plot');
   %  hold on
   % 
   %  % figure(7)
   %  % plot3(vx,vy,vz)
   %  % patch([vx' nan],[vy' nan],[vz' nan],[v_t1' nan],'EdgeColor','interp','FaceColor','none')
   %  % view(3)
   %  % title('3D plot of Velocity')
   %  % hold on
   %  % 
   %  % figure(8)
   %  % plot(vx,vz)
   %  % patch([vx' nan],[vz' nan],[v_t1' nan],'EdgeColor','interp','FaceColor','none')
   %  % title('V_X vs V_Z Plot')
   %  % hold on
   %  % 
   %  % figure(9)
   %  % plot(vy,vz)
   %  % patch([vy' nan],[vz' nan],[v_t1' nan],'EdgeColor','interp','FaceColor','none')
   %  % title('V_Y vs V_Z Plot')
   %  % hold on
   %  % 
   %  % figure(10)
   %  % plot(vx,vy)
   %  % patch([vx' nan],[vy' nan],[v_t1' nan],'EdgeColor','interp','FaceColor','none')
   %  % title('V_X vs V_Y Plot')
   %  % hold on
   %  % 
   %  figure(11)
   %  plot(v_t1,vz)
   %  patch([v_t1' nan],[vz' nan],[v_t1' nan],'EdgeColor','interp','FaceColor','none')
   %  title('V_Z vs Time Plot')
   %  hold on
   %  % 
   %  % figure(12)
   %  % plot(v_t1,vmag)
   %  % patch([v_t1' nan],[vmag' nan],[v_t1' nan],'EdgeColor','interp','FaceColor','none')
   %  % title('V_{mag} vs Time Plot')
   %  % hold on
   % 
   %  figure(13)
   %  plot3(x,y,z)
   %  patch([currTrace.row' nan],[currTrace.col' nan],[currTrace.z' nan],[t1' nan],'EdgeColor','interp','FaceColor','none')
   %  view(3)
   %  xlabel('X (\mum)','FontSize', 12, 'FontWeight', 'bold');
   %  ylabel('Y (\mum)','FontSize', 12, 'FontWeight', 'bold');
   %  zlabel('Z (\mum)','FontSize', 12, 'FontWeight', 'bold');
   %  xticks([0 20000 40000 60000])
   %  xticklabels({'0','20','40','60'})
   %  yticks([0 20000 40000 60000])
   %  yticklabels({'0','20','40','60'})
   %  %zticks([-2000 -1000 0 1000 2000])
   %  %zticklabels({'-2','-1','0','1','2'})
   %  set(gca,'fontweight','bold','LineWidth',1.5,'Box','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0]);
   %  ax = gca;
   %  ax.FontSize = 14; 
   %  ax.FontName='Times New Roman';
   %  %title('3D plot of position')
   %  hold on
   % 
   %  figure(14)
   %  plot(x,y)
   %  patch([x' nan],[y' nan],[t' nan],'EdgeColor','interp','FaceColor','none');
   %  %view(3)
   %  hold on
   % 
   %  figure(15)
   %  scatter(x(end),y(end),'filled');
   %  hold on
   % 
   %  figure(16)
   %  scatter(mean(x),mean(y),'filled');
   %  hold on

  %%%% End Plotting



end

x_std=mean(std_x)
y_std=mean(std_y)
z_std=mean(std_z)
z_mean=mean(mean_z)
range_x=max(max_x)-min(min_x);
range_y=max(max_y)-min(min_y);

%%%%%% Rolling std.

trap_height=zeros(length(GoodPos),1);
for i=1:length(GoodPos)

    currTrace=GoodPos{i,1};
   
    for j=1:height(currTrace)
        arr=currTrace.z(j:end);
        strd=std(arr);
        if(strd<=z_std/2  && length(arr)>=0.5*height(currTrace))
            trap_height(i)=mean(arr);
        end
    end
end
trap_height = trap_height(any(trap_height ~= 0, 2), :);
mean(trap_height)
std(trap_height)


%%%% && length(find(arr>=height_filter))>=0.50*length(arr)


%% cumulative trapping
check_n=0;
z_count=zeros(1,30);

for i=1:length(GoodPos)
    currTrace = GoodPos{i,1};
    x=currTrace.row;
    y=currTrace.col;
    z=currTrace.z;
    t=currTrace.rT/1000;

    % if test1 and test2 both true, then find the first point when
    % z>height_filter and x(and)y is in mask position
 
    test1=(length(find(z>height_filter)))>=0.50*height(currTrace); %half trace is above height filter
    
    n=0;
     for j=1:length(x)
         
         if(Mask1( round(x(j)), round(y(j)) )  == 1)
             n=n+1;
         end
         %check_n=[check_n; n/length(x)];  
         if(n>=0.8*length(x))             %% 70% of the trace is inside the mask
             test2=1;  % xy trace is in the masked zones
         else
             test=0;
         end
     end

  if(test1 && test2)
      p=find(z>= ((mean(trap_height))),1);
      if(~isempty(p))
      t1=t(p);
      t1=floor(t1);
      z_count(t1+1)=z_count(t1+1)+1;
      end
  end

end

z_count=z_count./num_videos;

%save('cumu_trap.mat','z_count');

figure
plot(cumsum(z_count),'o-');

figure
plot(z_count);

%% Z count

[z_count, edges] = histcounts(AllPos(:,3), 100);
binCenters = edges(1:end-1) + diff(edges)/2;

[a,b]=max(z_count);
binCenters(b)

figure
plot(binCenters, z_count,'o-', 'LineWidth', 2)
xlabel('Z (nm)')
ylabel('Count')

%% ideas

%better the mask %%%%
%story 
%release last 10s speed?? 
%diffusion compare speed with release