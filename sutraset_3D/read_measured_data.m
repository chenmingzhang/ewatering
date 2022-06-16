% Tbl = readtable('C:\Project\MDBA\data_deliverable\out.csv');
Tbl = readtable('/home/Student/s4680073/macondo/test/sp_sch.csv');

time                   = datetime(Tbl.date_time,'InputFormat','yyyy-MM-dd HH:mm:ss');
time_interval          = 600;
time_second            = (0:time_interval:(length(time)-1)*time_interval);
time_measured_data_day = time_second/86400;
% plot(time,Tbl.sa2_groundwater_table_mAHD)

% index for p c and s
x_idx  = strcmp(nod(1).label,'X');
y_idx  = strcmp(nod(1).label,'Y');
p_idx  = strcmp(nod(1).label,'Pressure');
c_idx  = strcmp(nod(1).label,'Concentration');
s_idx  = strcmp(nod(1).label,'Saturation');
sm_idx = strcmp(nod(1).label,'SM');
vy_idx = strcmp(ele(1).label,'Y velocity');
vx_idx = strcmp(ele(1).label,'X velocity');
rek_idx = strcmp(ele(1).label,'REK');
yele_idx = strcmp(ele(1).label,'Y origin');
%Now enter the path of model output
X          = zeros(length(nod),length(nod(1).terms{1,2}));
Y          = X;
Saturation = X;
t          = zeros(length(nod),1);

%for i=1:length(nod)
%    X(i,:)          = nod(i).terms{x_idx};
%    Y(i,:)          = nod(i).terms{y_idx};
%    Saturation(i,:) = nod(i).terms{};
%    t(i)            = nod(i).tout;
%end

sx        = linspace(0,150,151);  % For griddata of .nod data
sy        = linspace(-10,0,11);
sy(end+1) = 1;

[X,Y]                = meshgrid(sx,sy);
watertable_sa2_model = zeros(length(nod),1);
watertable_sa1_model = zeros(length(nod),1);

%for i=1:length(nod)
%    Saturation2d = reshape(Saturation(i,:),[length(sy),length(sx)]);
%    sa2_profile  = Saturation2d(:,1);
%    for j=1:length(sy)-1
%        if sa2_profile(j)==1 & sa2_profile(j+1)<1
%            watertable_sa2_model(i)=sy(j+1);
%        end
%    end
%    sa1_profile=Saturation2d(:,point3.x);
%        for j=1:length(sy)-1
%        if sa1_profile(j)==1 & sa1_profile(j+1)<1
%            watertable_sa1_model(i)=sy(j+1);
%        end
%    end
%%     contourf(X,Y,Saturation2d);
%%     pause; 
%end 
%watertable_diff_sa2_m=(watertable_sa2_model-watertable_sa2_model(1));
%watertable_diff_sa1_m=(watertable_sa1_model-watertable_sa1_model(1));


% plot(t(1:end),watertable_diff_m);
% hold on;
% plot(time_second,Tbl.sa2_groundwater_table_rise_m/1000)
