%% Pre-processing results - These values need to be changed, accordingly
start_time=2;
last_time=length(nod)-1;    % To see last time step
inp.nn1=105;                % Number of rows (.INP was too large to load)
inp.nn2=4001;               % Number of columns (.INP was too large to load)
inp.qinc=0.05;              % Thickness of surface cells 
sec=3600;
day=24;
conc_limit=0.0175;          % Concentration where XL and Vf are based upon
% Modified Analytical solution for Base case parameters proposed by Werner (2017): https://doi.org/10.1016/j.advwatres.2017.02.001
X=[200,	190,	180,	170,	160,	150,	140,	130,	120, 110, 100, 90, 80, 70];
Y=[-4.566544348	-3.56239023 -3.013656664 -2.586621068 -2.224529609 -1.904507484	-1.614603134 -1.347637756 -1.098899132 -0.865096615	-0.64382461	-0.433261333 -0.231987828 -0.038873304];

%% Create shortcuts
% Make x and y matrix (nodes)
x_matrix=reshape(nod(1).terms{xnod_idx},[inp.nn1,inp.nn2]);
y_matrix=reshape(nod(1).terms{ynod_idx},[inp.nn1,inp.nn2]);

% Make x and y matrix (elements)
xele_matrix=reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
yele_matrix=reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);

% Make pressure and sarutration array at top
p_top=arrayfun(@(y) y.terms{p_idx}(inp.nn1),nod);
sw_top=arrayfun(@(y) y.terms{s_idx}(inp.nn1),nod);

% Make Pressure matrix
% p_matrix = reshape(nod(nt).terms{p_idx},[inp.nn1,inp.nn2]); % Basic function
p_matrix_start = reshape(nod(start_time).terms{p_idx},[inp.nn1,inp.nn2]); % Initial pressure
p_matrix_end = reshape(nod(last_time).terms{p_idx},[inp.nn1,inp.nn2]); % Last Timestep

% Make Concentration matrix
% c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]); % Basic function
c_matrix_start = reshape(nod(start_time).terms{c_idx},[inp.nn1,inp.nn2]); % Initial concentration
c_matrix_end = reshape(nod(last_time).terms{c_idx},[inp.nn1,inp.nn2]); % Last Timestep

% Make Saturation matrix
%s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]); % Basic function
s_matrix_start  = reshape(nod(start_time).terms{s_idx},[inp.nn1,inp.nn2]); % Intial Saturation
s_matrix_end  = reshape(nod(last_time).terms{s_idx},[inp.nn1,inp.nn2]); % Last Timestep

% Make Velocity (x-direction) matrix (No initial velocity field)
%vx_matrix = reshape(ele(nt).terms{vx_idx}*24*3600,[inp.nn1-1,inp.nn2-1]); %Basic function (mm/d)
vx_matrix_start = reshape(ele(start_time).terms{vx_idx}*24*3600,[inp.nn1-1,inp.nn2-1]); % Timestep 1 (almost no velocity)
vx_matrix_end = reshape(ele(last_time-1).terms{vx_idx}*24*3600,[inp.nn1-1,inp.nn2-1]); % Last Timestep

% Make Velocity (y-direction) matrix (No initial velocity field)
%vy_matrix = reshape(ele(nt).terms{vy_idx}*24*3600,[inp.nn1-1,inp.nn2-1]); %Basic function (mm/d)
vy_matrix_start = reshape(ele(start_time).terms{vy_idx}*24*3600,[inp.nn1-1,inp.nn2-1]); % Timestep 1 (almost no velocity)
vy_matrix_end = reshape(ele(last_time-1).terms{vy_idx}*24*3600,[inp.nn1-1,inp.nn2-1]); % Last Timestep

% Create Area (Only for homogeneous cell field and where depth of cells is 1m)
Area            = zeros(inp.nn1,inp.nn2);
Area            = Area+(inp.qinc*inp.qinc);
Area(1,:)       = (inp.qinc/2)*inp.qinc;
Area(inp.nn1,:) = (inp.qinc/2)*inp.qinc;
Area(:,1)       = (inp.qinc/2)*inp.qinc;
Area(:,inp.nn2) = (inp.qinc/2)*inp.qinc;
Area(1,1)       = (inp.qinc/2)*(inp.qinc/2);Area(inp.nn1,1) = (inp.qinc/2)*(inp.qinc/2);Area(1,inp.nn2) = (inp.qinc/2)*(inp.qinc/2);Area(inp.nn1,inp.nn2) = (inp.qinc/2)*(inp.qinc/2);

%% Get INFORMATION
%Calculate Vf (Volume of freshwater with concentration smaller than 17.5 ppt)
c=reshape(nod(last_time).terms{c_idx},[inp.nn1,inp.nn2]);
Vf=0;
for j=1:inp.nn1
    for i=1:inp.nn2
        if c(j,i)<conc_limit
            Vf=Area(j,i)+Vf;
        end
    end
end
for j=1:length(bcop)
    for k=1:inp.nn2
        if c(1,k)<conc_limit
            fresh_top(k)=x_matrix(1,k);
        else
            fresh_top(k)=x_matrix(1,inp.nn2);
        end
        if c(inp.nn1,k)<conc_limit
            fresh_bottom(k)=x_matrix(inp.nn1,k);
        else
            fresh_bottom(k)=x_matrix(1,inp.nn2);
        end
    end
    fresh_extent_top(j)=x_matrix(1,inp.nn2)-min(fresh_top);
    fresh_extent_bottom(j)=x_matrix(1,inp.nn2)-min(fresh_bottom);
    time(j)=bcop(j).tout/3600/24/365;
end
fresh_extent_top=fresh_extent_top.';
fresh_extent_bottom=fresh_extent_bottom.';
time=time.';
tab2=table(time,fresh_extent_top,fresh_extent_bottom);
na=['info.txt'];
writetable(tab2,na,'Delimiter',' '); 

% Waterbalance
for i=length(bcop);  
    left_ef=0;
    left_in=0;
    riv_ef=0;
    riv_in=0;
    for j=1:inp.nn1
        if bcop(i).qpl(j)<0
            riv_ef=riv_ef+bcop(i).qpl(j);
        else
            riv_in=riv_in+bcop(i).qpl(j);
        end
        if bcop(i).qpl(j+inp.nn1)<0
            left_ef=left_ef+bcop(i).qpl(j+inp.nn1);
        else
            left_in=left_in+bcop(i).qpl(j+inp.nn1);
        end
    end
end
Left_w_efflux=left_ef*-day*sec;
Left_w_influx=left_in*day*sec;
River_w_efflux=riv_ef*-day*sec;
River_w_influx=riv_in*day*sec;
Total_w_influx=Left_w_influx+River_w_influx;
Total_w_efflux=Left_w_efflux+River_w_efflux;


p_Left_w_efflux=((Total_w_efflux-Left_w_efflux)/Total_w_efflux)*100;
p_Left_w_influx=((Total_w_influx-Left_w_influx)/Total_w_influx)*100;
p_River_w_efflux=((Total_w_efflux-River_w_efflux)/Total_w_efflux)*100;
p_River_w_influx=((Total_w_influx-River_w_influx)/Total_w_influx)*100;
clear left_ef left_in riv_ef riv_in

for i=length(bcop);  
    left_ef=0;
    left_in=0;
    riv_ef=0;
    riv_in=0;
    for j=1:inp.nn1
        if bcop(i).qpu(j)<0
            riv_ef=riv_ef+bcop(i).qpu(j);
        else
            riv_in=riv_in+bcop(i).qpu(j);
        end
        if bcop(i).qpu(j+inp.nn1)<0
            left_ef=left_ef+bcop(i).qpu(j+inp.nn1);
        else
            left_in=left_in+bcop(i).qpu(j+inp.nn1);
        end
    end
end
Left_s_efflux=left_ef*-day*sec;
Left_s_influx=left_in*day*sec;
River_s_efflux=riv_ef*-day*sec;
River_s_influx=riv_in*day*sec;
Total_s_influx=Left_s_influx+River_s_influx;
Total_s_efflux=Left_s_efflux+River_s_efflux;


p_Left_s_efflux=((Total_s_efflux-Left_s_efflux)/Total_s_efflux)*100;
p_Left_s_influx=((Total_s_influx-Left_s_influx)/Total_s_influx)*100;
p_River_s_efflux=((Total_s_efflux-River_s_efflux)/Total_s_efflux)*100;
p_River_s_influx=((Total_s_influx-River_s_influx)/Total_s_influx)*100;
clear left_ef left_in riv_ef riv_in

%Make table
Waterbalance = {'River';'Left boundary';'Total'};
Water_Influx = [River_w_influx;Left_w_influx;Total_w_influx];
Water_Efflux = [River_w_efflux;Left_w_efflux;Total_w_efflux];  
Salt_Influx = [River_s_influx;Left_s_influx;Total_s_influx];
Salt_Efflux = [River_s_efflux;Left_s_efflux;Total_s_efflux];
T_waterbalance = table(Waterbalance,Water_Influx,Water_Efflux,Salt_Influx,Salt_Efflux); 
T_waterbalance;
na=['waterbalance.txt'];
writetable(T_waterbalance,na,'Delimiter',' ');  
%Make table-Percentage
Waterbalance = {'River';'Left boundary'};
Water_Influx = [p_River_w_influx;p_Left_w_influx];
Water_Efflux = [p_River_w_efflux;p_Left_w_efflux];  
Salt_Influx = [p_River_s_influx;p_Left_s_influx];
Salt_Efflux = [p_River_s_efflux;p_Left_s_efflux];
T_waterbalance = table(Waterbalance,Water_Influx,Water_Efflux,Salt_Influx,Salt_Efflux); 
T_waterbalance;
na=['Percentage_waterbalance.txt'];
writetable(T_waterbalance,na,'Delimiter',' '); 

%% Creat Figures
% Plot pressure, concentration with velocity field and saturation contours in 1 figure

figure
ax1=subplot(3,2,1)              %Pressure timestep 1
contourf(x_matrix,y_matrix,p_matrix_start);
colormap(ax1,parula)
title('Pressure, timestep 1')
xlabel('Distance x (m)')
ylabel('Depth (m)')
b=colorbar;
b.Label.String ='Pressure (Pa)';

ax2=subplot(3,2,2)              %Pressure last timestep
contourf(x_matrix,y_matrix,p_matrix_end);
colormap(ax2,parula)
title('Pressure, last timestep')
xlabel('Distance x (m)')
ylabel('Depth (m)')
b=colorbar;
b.Label.String ='Pressure (Pa)';

ax3=subplot(3,2,3) ;            %Concentration timestep 1
contourf(x_matrix,y_matrix,c_matrix_start);
colormap(ax3,jet)
title('Concentration, timestep 1')
xlabel('Distance x (m)')
ylabel('Depth (m)')
a=colorbar;
a.Label.String ='Concentration (-)';
hold on
%quiver(xele_matrix,yele_matrix,vx_matrix_start,vy_matrix_start, 'color',[1 1 1],'AutoScaleFactor',0.5)

ax4=subplot(3,2,4) ;            %Concentration last timestep
contourf(x_matrix,y_matrix,c_matrix_end);
colormap(ax4,jet)
title('Concentration, last timestep')
xlabel('Distance x (m)')
ylabel('Depth (m)')
a=colorbar;
a.Label.String ='Concentration (-)';
hold on
%quiver(xele_matrix,yele_matrix,vx_matrix_end,vy_matrix_end, 'color',[1 1 1],'AutoScaleFactor',0.5)

ax5=subplot(3,2,5);             %Saturation timestep 1    
contourf(x_matrix(1:3,:),y_matrix(1:3,:),s_matrix_start(1:3,:)); 
colormap(ax5,parula);
title('Saturation, timestep 1')
xlabel('Distance x (m)')
ylabel('Depth (m)')
a=colorbar;
a.Label.String ='Saturation (-)';

ax6=subplot(3,2,6);             %Saturation last timestep  
contourf(x_matrix(1:3,:),y_matrix(1:3,:),s_matrix_end(1:3,:)); 
colormap(ax6,parula);
title('Saturation, last timestep')
xlabel('Distance x (m)')
ylabel('Depth (m)')
a=colorbar;
a.Label.String ='Saturation (-)';
hold on
print(gcf,['Overview.png'],'-dpng','-r300');

%Concentration last timestep
figure
ax2=subplot('Position',[0.08 0.7 0.84 0.2]);
contourf(x_matrix,y_matrix,log10(c_matrix_end),'LevelStep',0.005,'LineStyle','None');hold on
set(gca,'ylim',[-5.2 0.0],'YTick',[-5.2 -4.2 -3.2 -2.2 -1.2 -0.2],'YTicklabel',[0 1 2 3 4 5],'xlim',[0 200],'XTick',[0 20 40 60 80 100 120 140 160 180 200],'XTicklabel',[200 180 160 140 120 100 80 60 40 20 0],'fontsize',10)
xlabel('x (m)')
ylabel('z (m)')
caxis([-4 -1.45]);
colormap(ax2,jet);
c = colorbar('location','south','position',[0.4 0.57 0.52 0.02]);
c.Ticks = [-4 -3 -1.45];
c.TickLabels={'0','1','35'};
c.Label.String ='Concentration (ppt)';
c.Label.FontSize = 8;
hold on
plot(X,Y,'k--','LineWidth',2)
print(gcf,['Concentration.png'],'-dpng','-r300');
