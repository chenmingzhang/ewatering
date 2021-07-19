%et1_kgs   = -arrayfun(@(y) y.qin(1),bcof);
%area1     = 0.5D-6; % hard coding
%et1_mmday = et1_kgs/area1*86400;
%time_day  = [bcof.tout]/3600/24;
time_nod_day= arrayfun(@(y) y.tout,nod) * c.dayPsec;


% saturation profile
%d(1).terms{p_idx}(1:inp.nn1)

%dp0_idx=min(0-nod(1).terms{2}(1:inp.nn1)); 

% write x and y coordinates in matrix form.
x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);


dp0_idx = inp.nn1;
dp1_idx = inp.nn1-5;

p_top  = arrayfun(@(y) y.terms{p_idx}(inp.nn1),nod);
sw_top = arrayfun(@(y) y.terms{s_idx}(inp.nn1),nod);


fig_pos.left   = 0.05;
fig_pos.bottom = 0.8;
fig_pos.length = 0.3;
fig_pos.height = 0.15;
fig_pos.xy_plot_length = 0.15;

%nt=10;
%a.fig=figure;
a.fs = 10;
a.lw = 2;
a.cz = 8;
fs   = 5; % sampling frequency
% Creating the movie : quality = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt = 100;
%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig = figure;
%set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]);  % maximize the plotting figure
set(gcf,'Units','normalized', 'OuterPosition',[0 0 1 1]);  % maximize the plotting figure
mov           =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
mov.FrameRate = 5;mov.Quality=qt;
open(mov);




point1.x           = 5;
point1.y           = -6;
xdist_array        = point1.x-x_nod_array;
ydist_array        = point1.y-y_nod_array;
dist_array         = (xdist_array.^2+ydist_array.^2).^0.5;
[~,point1.idx_nod] = min(dist_array);
point1.x_nod_m     = x_nod_array(point1.idx_nod);  % the actual location of the node
point1.y_nod_m     = y_nod_array(point1.idx_nod);
point1.mark_shape       = 'rx';
point1.legend      = sprintf('x = %i ; y = %i',point1.x,point1.y);

point2.x           = 100;
point2.y           = -6;
xdist_array        = point2.x-x_nod_array;
ydist_array        = point2.y-y_nod_array;
dist_array         = (xdist_array.^2+ydist_array.^2).^0.5;
[~,point2.idx_nod] = min(dist_array);
point2.x_nod_m     = x_nod_array(point2.idx_nod);
point2.y_nod_m     = y_nod_array(point2.idx_nod);
point2.mark_shape       = 'gx';
point2.legend      = sprintf('x = %i ; y = %i',point2.x,point2.y);

point3.x           = 120;
point3.y           = -6;
xdist_array        = point3.x-x_nod_array; ydist_array  = point3.y-y_nod_array;
dist_array         = (xdist_array.^2+ydist_array.^2).^0.5;
[~,point3.idx_nod] = min(dist_array);
point3.x_nod_m     = x_nod_array(point3.idx_nod);
point3.y_nod_m     = y_nod_array(point3.idx_nod);
point3.legend      = sprintf('x= %i ; y= %i',point3.x,point3.y);
point3.mark_shape       = 'bx';

% samping a pond
point_pond.x           = 5;
point_pond.y           = 0;
xdist_array            = point_pond.x-x_nod_array; ydist_array = point_pond.y-y_nod_array;
dist_array             = (xdist_array.^2+ydist_array.^2).^0.5;
[~,point_pond.idx_nod] = min(dist_array);
point_pond.x_nod_m     = x_nod_array(point_pond.idx_nod);
point_pond.y_nod_m     = y_nod_array(point_pond.idx_nod);
point_pond.legend      = sprintf('x = %i ; y = %i',point_pond.x,point_pond.y);
point_pond.mark_shape       = 'cx';


vertical_line_pond_centre.x       = 1.0;
abs_xdist_array_m                 = abs(vertical_line_pond_centre.x-x_nod_array);
min_abs_xdist                     = min(abs_xdist_array_m);
vertical_line_pond_centre.idx_nod = find( abs_xdist_array_m == min_abs_xdist);


vertical_line_pond_boundary.x       = pond_radius_m;    %50
abs_xdist_array_m                   = abs(vertical_line_pond_boundary.x-x_nod_array);
min_abs_xdist                       = min(abs_xdist_array_m);
vertical_line_pond_boundary.idx_nod = find( abs_xdist_array_m == min_abs_xdist);


vertical_line_pond_outside.x       = pond_radius_m+10;    % 60m
abs_xdist_array_m                   = abs(vertical_line_pond_outside.x-x_nod_array);
min_abs_xdist                       = min(abs_xdist_array_m);
vertical_line_pond_outside.idx_nod = find( abs_xdist_array_m == min_abs_xdist);

for nt=2:2:length(nod)-1

    fprintf('plotting the %i out of %i results,\n', nt, length(nod));
    number_of_negative_c=sum(nod(nt).terms{c_idx}<0);
    [max_c,max_c_idx]=max(nod(nt).terms{c_idx});
    fprintf('    max concentration is %f at location  %i, x= %f , y= %f \n',max_c,max_c_idx,  x_nod_mtx(max_c_idx),  y_nod_mtx(max_c_idx));
    if number_of_negative_c>0;
          fprintf('    %i negative concentration\n', length(number_of_negative_c));
    end

    %% -------------  sub 1 ET over time  --------------
    a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom,...
          fig_pos.length,fig_pos.height]);
    %a.plot1=plot(time_day(1:nt),et1_mmday(1:nt),...
    %         'k-','linewidth',a.lw);hold on
    %%a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    %get(gca,'xtick');
    %set(gca,'fontsize',10);
    %txt=sprintf('Result at day %.2f , k_silt = %.2e m2, k_sand = %.2e m2, init P in silt layer = %.2e Pa'...
    %    ,nod(nt).tout*c.dayPsec,permeability_silt_m2,permeability_sand_m2,initial_silt_layer_pressure_pa);
    %title(txt, 'Interpreter', 'none');

    %xlabel('Time (day)','FontSize',a.fs);
    %ylabel('Evt (mm/day)','FontSize',a.fs);
    %axis([-1 time_day(end) -1 21])
    %title('ead profile');
    %legend('show','Location','East')
    
    
    
    %% -------------  sub 2 sat over time  --------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+0.5,fig_pos.bottom,...
          fig_pos.length,fig_pos.height]);
    %arrayfun(@(y) y.tout,nod(1:nt))

    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{s_idx}(point1.idx_nod),nod(1:nt)),...
             'r-','linewidth',a.lw);hold on
    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{s_idx}(point1.idx_nod),nod(1:nt)),...
             'g-','linewidth',a.lw);hold on
    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{s_idx}(point3.idx_nod),nod(1:nt)),...
             'b-','linewidth',a.lw);hold on
    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{s_idx}(point_pond.idx_nod),nod(1:nt)),...
             'c-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel('Saturation (-)','FontSize',a.fs);
    axis([-1 time_nod_day(end) -0.1 1.1])
    
    
    %% -------------  sub 3 conc over time  --------------
    a.sub3=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{c_idx}(point1.idx_nod),nod(1:nt)),...
             'r-','linewidth',a.lw);hold on
    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{c_idx}(point2.idx_nod),nod(1:nt)),...
             'g-','linewidth',a.lw);hold on
    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{c_idx}(point3.idx_nod),nod(1:nt)),...
             'b-','linewidth',a.lw);hold on
    a.plot2=plot(time_nod_day(1:nt),...
           arrayfun(@(y) y.terms{c_idx}(point_pond.idx_nod),nod(1:nt)),...
             'c-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel('Salinity (kg/kg)','FontSize',a.fs);
    axis([-1 time_nod_day(end) -0.005 0.04])
    
    
    %% -------------  sub 4 pressure over time  --------------
    
    %salt_dep=arrayfun(@(y) y.terms{sm_idx}(inp.nn1),nod(1:nt))/2165/area1*1000;
    a.sub4=subplot('position'...
         ,[fig_pos.left+0.5,fig_pos.bottom-fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    a.plot2=plot(time_nod_day(1:nt),...
            arrayfun(@(y) y.terms{p_idx}(point1.idx_nod),nod(1:nt)), ...  
             'r-','linewidth',a.lw,'displayname',point1.legend);hold on
    a.plot2=plot(time_nod_day(1:nt),...
            arrayfun(@(y) y.terms{p_idx}(point2.idx_nod),nod(1:nt)), ...  
             'g-','linewidth',a.lw,'displayname',point2.legend);hold on
    a.plot2=plot(time_nod_day(1:nt),...
            arrayfun(@(y) y.terms{p_idx}(point3.idx_nod),nod(1:nt)), ...  
             'b-','linewidth',a.lw,'displayname',point3.legend);hold on
    a.plot2=plot(time_nod_day(1:nt),...
            arrayfun(@(y) y.terms{p_idx}(point_pond.idx_nod),nod(1:nt)), ...  
             'c-','linewidth',a.lw,'displayname',point3.legend);hold off
    get(gca,'xtick');
    set(gca,'fontsize',10);
    legend('location','east','fontsize',5)
    xlabel('Time (day)','FontSize',a.fs);
    ylabel('pressure (pa)','FontSize',a.fs);
    %axis([-1 time_day(end) -0.005 0.3])
    
    
    %% -------- contour plot on Saturation ---------
    a.sub5=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-2*fig_pos.height-0.07,...
          fig_pos.length+0.02,fig_pos.height]);
    % write pressure and conc in matrix form.
    s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);
    
    contourf(x_matrix,y_matrix,s_matrix);hold on
    colormap(a.sub5,parula);
    scatter(point1.x,point1.y,point1.mark_shape,'linewidth',2);hold on
    scatter(point2.x,point2.y,point2.mark_shape,'linewidth',2);hold on
    scatter(point3.x,point3.y,point3.mark_shape,'linewidth',2);hold on
    scatter(point_pond.x,point_pond.y,point_pond.mark_shape,'linewidth',2);hold off
    colorbar
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Saturation')
    
    
    %% -------- contour plot on concentration ---------
    a.sub6=subplot('position'...
         ,[fig_pos.left+0.5,...
         fig_pos.bottom-2*fig_pos.height-0.07,...
          fig_pos.length+0.02,fig_pos.height]);
    
    % write pressure and conc in matrix form.
    c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]);
    
    contourf(x_matrix,y_matrix,c_matrix);hold on
    colormap(a.sub6,jet);
    scatter(point1.x,point1.y,point1.mark_shape,'linewidth',2);hold on
    scatter(point2.x,point2.y,point2.mark_shape,'linewidth',2);hold on
    scatter(point3.x,point3.y,point3.mark_shape,'linewidth',2);hold on
    scatter(point_pond.x,point_pond.y,point_pond.mark_shape,'linewidth',2);hold off
    colorbar
    get(gca,'xtick');
    set(gca,'fontsize',10);
    title('Salinity (kg/kg)');
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    %axis([10, 40,9,10])
    
    
    %% --------- plot saturation profiles  ------------------
    a.sub7=subplot('position'...
         ,[fig_pos.left,...
         fig_pos.bottom-3*fig_pos.height-0.25,...
          fig_pos.xy_plot_length,fig_pos.height+0.1]);
    plot(nod(nt).terms{s_idx}(vertical_line_pond_centre.idx_nod)...
         ,nod(nt).terms{y_idx}(vertical_line_pond_centre.idx_nod)...
         ,'r-','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{s_idx}(vertical_line_pond_boundary.idx_nod)...
         ,nod(nt).terms{y_idx}(vertical_line_pond_boundary.idx_nod)...
         ,'g-','linewidth',a.lw) ;hold on

    plot(nod(nt).terms{s_idx}(vertical_line_pond_outside.idx_nod)...
         ,nod(nt).terms{y_idx}(vertical_line_pond_outside.idx_nod)...
         ,'b-','linewidth',a.lw) ;hold on

    %plot(slab(2,:,1),slab(1,:,1),'rd',slab(2,:,2),...
    %    slab(1,:,2),'go',slab(2,:,3),slab(1,:,3),'cx','linewidth',a.lw);hold off;
    %axis([-0.05, 1.05,0,0.05])
    x_rectangle = [0,1,1,0];
    y_rectangle = [-6,-6,-10,-10];
    p=patch(x_rectangle,y_rectangle,'r');hold off
    set(p,'FaceAlpha',0.5);
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('Saturation (-)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    
    %% --------- plot conc profiles  ------------------
    a.sub8=subplot('position'...
         ,[fig_pos.left+(fig_pos.xy_plot_length+0.05),...
         fig_pos.bottom-3*fig_pos.height-0.25,...
          fig_pos.xy_plot_length,fig_pos.height+0.1]);

    plot(nod(nt).terms{c_idx}(vertical_line_pond_centre.idx_nod)...
         ,nod(nt).terms{y_idx}(vertical_line_pond_centre.idx_nod)...
         ,'r-','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{c_idx}(vertical_line_pond_boundary.idx_nod)...
         ,nod(nt).terms{y_idx}(vertical_line_pond_boundary.idx_nod)...
         ,'g-','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{c_idx}(vertical_line_pond_outside.idx_nod)...
         ,nod(nt).terms{y_idx}(vertical_line_pond_outside.idx_nod)...
         ,'b-','linewidth',a.lw) ;hold on
    x_rectangle = [0,0.04,0.04,0];
    y_rectangle = [-6,-6,-10,-10];
    p=patch(x_rectangle,y_rectangle,'r');hold off
    set(p,'FaceAlpha',0.5);
    %rectangle('Position',[0 -5 1 5],'FaceColor','red','FaceAlpha',0.5);  % no facealpha
    %plot(clab(2,:,1),clab(1,:,1),'rd',clab(2,:,2),clab(1,:,2)...
    %,'go',clab(2,:,3),clab(1,:,3),'cx','MarkerSize',a.cz,'linewidth',a.lw);hold off
    %axis([-0.05, 0.28,0,0.05])
    get(gca,'xtick');
    set(gca,'fontsize',12);
    xlabel('Salinity (kg/kg)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    
    %% --------- plot x velocity  ------------------
    a.sub9=subplot('position'...
         ,[fig_pos.left+2*(fig_pos.xy_plot_length+0.05),...
         fig_pos.bottom-3*fig_pos.height-0.25,...
          fig_pos.xy_plot_length,fig_pos.height+0.1]);
    plot(ele(nt).terms{vx_idx}(1:inp.nn1-1)*c.secPday ...
         ,ele(nt).terms{yele_idx}(1:inp.nn1-1)...
         ,'linewidth',a.lw) ;hold off
    
    %% --------- plot relative permeability velocity  ------------------
    a.sub9=subplot('position'...
         ,[fig_pos.left+3*(fig_pos.xy_plot_length+0.05),...
         fig_pos.bottom-3*fig_pos.height-0.25,...
          fig_pos.xy_plot_length,fig_pos.height+0.1]);

    semilogx(ele(nt).terms{rek_idx}(1:inp.nn1-1) ...
         ,ele(nt).terms{yele_idx}(1:inp.nn1-1)...
         ,'linewidth',a.lw) ;hold off

    get(gca,'xtick');
    set(gca,'fontsize',12);
    xlabel('Relative permeability (-)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    
    saveas(gcf,['chart',num2str(nt),'.png'])
    F = getframe(gcf); % save the current figure
    writeVideo(mov,F);% add it as the next frame of the movie
end

close(mov);
close(a.fig);
