itout                                                            = extractfield(bcof,'itout');
durn                                                             = extractfield(bcof,'durn');
tout_s                                                           = extractfield(bcof,'tout');

dx_ele_mtx_3d                                                    = cat(3,zeros(nez,ney,1),diff(x_ele_mtx_3d,1,3) );
dz_ele_mtx_3d                                                    = cat(1,zeros(1,ney,nex),diff(z_ele_mtx_3d,1,1) );

area_xz_cell_mtx_m2_3d                                           = -dx_cell_mtx_3d .* dz_cell_mtx_3d;
area_xz_ele_mtx_m2_3d                                            = -dx_ele_mtx_3d  .* dz_ele_mtx_3d ;
area_xy_ele_mtx_m2_2d(1,:)                                       = (1/2) * (y_nod_matrix_mtx_2d(1,1:end-1)+y_nod_matrix_mtx_2d(1,2:end) ) .* squeeze(dx_cell_mtx_3d(1,3,2:end) )';
area_xy_ele_mtx_m2_2d(2,:)                                       = area_xy_ele_mtx_m2_2d(1,:);

for i = 1:nez
    area_xy_ele_mtx_m2_3d(i,:,:)                                 = area_xy_ele_mtx_m2_2d;
end
% source_sink_kgPs                                               =extractfield(bcof,'qin');
% source_sink_kgPs                                               =source_sink_kgPs(~isnan(source_sink_kgPs) );
% source_sink_kgPs_mtx                                           =reshape(source_sink_kgPs,[length(inp.iqcp) length(bcof)]);
qin_kgPs                                                         = cell2mat(arrayfun(@(y) y.qin(1:length(source_point_index) ),bcof,'UniformOutput',false) );
sink_kgPs                                                        = cell2mat(arrayfun(@(y) y.qin(length(source_point_index) + 1:end),bcof,'UniformOutput',false) );
qpl                                                              = cell2mat( arrayfun( @(y) y.qpl,bcop,'UniformOutput',false ) );
volume_cell_mtx_m3                                               = dx_cell_mtx_3d .* abs(dy_cell_mtx_3d) .* dz_cell_mtx_3d;
tout_diff_s                                                      = diff(tout_s);
tout_diff_s                                                      = [0 tout_diff_s];
input_diff_kg                                                    = ones(length(source_point_index),1) * tout_diff_s .* qin_kgPs;
input_total_kg                                                   = cumsum(sum(input_diff_kg,1) );

outflow_bcop_node_kg                                             = qpl .* ([itout(1) diff(itout)] .* durn);
outflow_bcop_node_sum_kg                                         = sum(outflow_bcop_node_kg,1);
%      outflow_bcop_sum_kg(i)=sum(c.rhow_pure_water*outflow_bcop_node_time_mtx_kg(:,i)./(c.rhow_pure_water+inp.drwdu*bcop(i).uucut) );
outflow_bcop_cumsum_kg                                           = cumsum(outflow_bcop_node_sum_kg);
% qpl                                                            =cell2mat(arrayfun(@(y) y.qpl,bcop,'UniformOutput',false) );
%et1_kgs                                                         = -arrayfun(@(y) y.qin(1),bcof);
%area1                                                           = 0.5D-6; % hard coding
%et1_mmday                                                       = et1_kgs/area1*86400;
%time_day                                                        = [bcof.tout]/3600/24;
time_nod_day                                                     = arrayfun(@(y) y.tout,nod) * c.dayPsec;
time_series                                                      =([itout(1) diff(itout)] .* durn);
pet_observed_kg                                                  = cumsum(time_interval*c.rhow_pure_water*Tbl.pet_mmPday*c.mPmm*(pi*pond_radius_m^2)*c.dayPsec);
% saturation profile
%d(1).terms{p_idx}(1:inp.nn1)
%calculate evaporated water from the sink point
et_mmPday                                                        =c.m2mm*sink_kgPs/c.rhow_pure_water./evap_point_area_xy_array/c.dayPsec;
%%calculate movmean of pet from observation data
pet_movmean_mmPday                                               =movmean(Tbl.pet_mmPday,6*24);
%calculate flux cross the groundwater table
saturation                                                       = cell2mat(arrayfun(@(y) y.terms{s_idx},nod,'UniformOutput',false) );
vy                                                               = cell2mat(arrayfun(@(y) y.terms{vy_idx},ele,'UniformOutput',false) );
vy_matrix                                                        = reshape(vy,[inp.nn1-1,inp.nn2-1,inp.nn3-1,length(ele)]);
vz                                                               = cell2mat(arrayfun(@(y) y.terms{vz_idx},ele,'UniformOutput',false) );
vz_matrix                                                        = reshape(vz,[inp.nn1-1,inp.nn2-1,inp.nn3-1,length(ele)]);
saturation_xyzt                                                  = reshape(saturation,[inp.nn1,inp.nn2,inp.nn3,length(nod)]);

pond_boundary_ele_index                                          = find(mask_ele_mtx_pond_boundary == 1) ;
mask_ele_mtx_pond_boundary_for_flux                              = zeros(nez,ney,nex);
mask_ele_mtx_pond_boundary_for_flux(pond_boundary_ele_index+2)   = 1;
mask_ele_mtx_pond_boundary_for_flux                              = logical(mask_ele_mtx_pond_boundary_for_flux);
vz_matrix_for_flux                                               = vz_matrix;
vz_matrix_for_flux(vz_matrix_for_flux>0)                         = 0;
for i =1:length(ele)
    flux_ele_xy_kgPs(:,:,:,i)                                    = area_xy_ele_mtx_m2_3d .* nod_to_ele(nex,ney,nez,squeeze(saturation_xyzt(:,:,:,i) ) )...
																	 .* nod_to_ele(nex,ney,nez,porosity_nod_mtx) .* -squeeze(vz_matrix_for_flux(:,:,:,i) ) .* (c.rhow_pure_water);%flux across xy plane
    flux_ele_xz_kgPs(:,:,:,i)                                    = area_xz_ele_mtx_m2_3d .* nod_to_ele(nex,ney,nez,squeeze(saturation_xyzt(:,:,:,i) ) )...
																 .* nod_to_ele(nex,ney,nez,porosity_nod_mtx) .* -squeeze(vy_matrix(:,:,:,i) ) .* (c.rhow_pure_water);%flux across xy plane       
    pond_flux_ele_xy_kgPs(:,:,:,i)                               = squeeze(flux_ele_xy_kgPs(:,:,:,i) ) .* mask_ele_mtx_pond_boundary_for_flux;
    aquitard_flux_ele_xy_kgPs(:,:,:,i)                           = squeeze(flux_ele_xy_kgPs(:,:,:,i) ) .* mask_ele_mtx_aquitard_layer;%%only consider the downward flow
    silt_flux_ele_xz_kgPs(:,:,:,i)                               = squeeze(flux_ele_xz_kgPs(:,:,:,i) ) .* (mask_ele_mtx_silt_layer - mask_ele_mtx_pond_layer);
end

Infiltration_kg                                                  = cumsum( squeeze(sum(sum(sum( (pond_flux_ele_xy_kgPs),1) , 2 ) , 3) ) )' .* tout_diff_s;
Infiltration_macropore_kg                                        = cumsum( squeeze(sum(sum (squeeze(pond_flux_ele_xy_kgPs(:,1,:,:) ),1) , 2 ) ) )' .* tout_diff_s;
Infiltration_matrix_kg                                           = cumsum( squeeze(sum(sum (squeeze(pond_flux_ele_xy_kgPs(:,2,:,:) ),1) , 2 ) ) )' .* tout_diff_s;

Recharge_kg                                                      = cumsum( squeeze(sum(sum(sum( (aquitard_flux_ele_xy_kgPs),1) , 2 ) , 3) ) )' .* tout_diff_s;
Recharge_macropore_kg                                            = cumsum( squeeze(sum(sum (squeeze(aquitard_flux_ele_xy_kgPs(:,1,:,:) ),1) , 2 ) ) )' .* tout_diff_s;
Recharge_matrix_kg                                               = cumsum( squeeze(sum(sum (squeeze(aquitard_flux_ele_xy_kgPs(:,2,:,:) ),1) , 2 ) ) )' .* tout_diff_s;
exchange_macropore_kg                                            = -cumsum( squeeze(sum(sum( squeeze(silt_flux_ele_xz_kgPs(:,1,:,:) ),1 ), 2 ) ) )' .* tout_diff_s;
exchange_matrix_kg                                               = -cumsum( squeeze(sum(sum( squeeze(silt_flux_ele_xz_kgPs(:,2,:,:) ),1 ), 2 ) ) )' .* tout_diff_s;

% pond_flux_xy_kgPs(i)                                           = squeeze(-sum(sum(sum(pond_flux_ele_xy_kgPs,1),2),3) );
for i =1:length(nod)
    s_matrix                                                     = reshape(nod(i).terms{s_idx},[inp.nn1,inp.nn2,inp.nn3]);   
    water_mass_per_cell_kg                                       = volume_cell_mtx_m3 .* porosity_nod_mtx .* s_matrix * ( c.rhow_pure_water );
    water_height_per_cell_mtx_m                                  = dy_cell_mtx_3d .* porosity_nod_mtx .* s_matrix;
    % water_amount_kg                                            =volume_cell_m3 .* inp.por .* saturation .* (c.rhow_pure_water);
    % pond_water_amount_kg                                       =volume_cell_m3 .* porosity_nod_mtx_gravity_compensated_1d .* saturation_gravity_compensated_1d .* (c.rhow_pure_water+inp.drwdu*solute_concentration_gravity_compensated_1d) .* reshape(mask_nod_mtx_pond_cell_gravity_compensated,[inp.nn1*inp.nn2,1]);
    pond_water_cell_kg                                           = water_mass_per_cell_kg .* mask_nod_mtx_pond_cell;
    above_groundwater_zone_water_cell_kg                         = water_mass_per_cell_kg .* (mask_porosity_nod_silt_layer-mask_nod_mtx_macropore-mask_nod_mtx_pond_cell);
    total_water_sum_kg(i)                                        = squeeze(sum(sum(sum(water_mass_per_cell_kg,1 ),2),3) )';
    pond_water_sum_kg(i)                                         = squeeze(sum(sum(sum(pond_water_cell_kg,1),2),3) )';
    above_groundwater_zone_water_cell_kg                         = water_mass_per_cell_kg .* (mask_nod_mtx_silt_layer - mask_nod_mtx_pond_cell);
    silt_zone_water_cell_kg(i)                                   = squeeze(sum(sum(sum(above_groundwater_zone_water_cell_kg,1),2),3 ) );
end
% write x y and z coordinates in matrix form.
water_height_per_cell_mtx_m                                      = -dz_cell_mtx_3d .* saturation_xyzt;

porosity_ele_mtx                                                 = nod_to_ele(nex,ney,nez,porosity_nod_mtx);
% p_top                                                          = arrayfun(@(y) y.terms{p_idx}(inp.nn1),nod);
% sw_top                                                         = arrayfun(@(y) y.terms{s_idx}(inp.nn1),nod);
point_1_pond.x                                                   = 5;
[M index_point_1_pond]                                           = min(abs(point_1_pond.x-x_array) );
point_1_pond.y                                                   = 2 * pi * point_1_pond.x;
point_1_pond.z                                                   = location_of_pond_boundary_nod(index_point_1_pond);
xdist_array                                                      = point_1_pond.x-x_nod_array; ydist_array = point_1_pond.y-y_nod_array;zdist_array = point_1_pond.z - z_nod_array;
dist_array                                                       = (xdist_array .^ 2+ydist_array .^ 2+zdist_array .^ 2) .^ 0.5;
[~,point_1_pond.idx_nod]                                         = min(dist_array);
point_1_pond.x_nod_m                                             = x_nod_array(point_1_pond.idx_nod);
point_1_pond.y_nod_m                                             = y_nod_array(point_1_pond.idx_nod);
point_1_pond.z_nod_m                                             = z_nod_array(point_1_pond.idx_nod);
point_1_pond.legend                                              = sprintf('x = %i ; y = %i; z = %i',point_1_pond.x,point_1_pond.z);
point_1_pond.mark_shape                                          = 'cx';

%%
sa2_macropore_point_a.x                                          = 5;
[M index_sa2_macropore_point_a]                                  = min(abs(sa2_macropore_point_a.x-x_array) );
sa2_macropore_point_a.y                                          = -2 * pi * sa2_macropore_point_a.x;
sa2_macropore_point_a.z                                          = depth_groundwater_table_m + 1;
xdist_array                                                      = sa2_macropore_point_a.x-x_nod_array; ydist_array = sa2_macropore_point_a.y-y_nod_array;zdist_array = sa2_macropore_point_a.z-z_nod_array;
dist_array                                                       = (xdist_array .^ 2+ydist_array  .^ 2 +zdist_array  .^ 2 ).^0.5;
[~,sa2_macropore_point_a.idx_nod]                                = min(dist_array);
sa2_macropore_point_a.x_nod_m                                    = x_nod_array(sa2_macropore_point_a.idx_nod);
sa2_macropore_point_a.y_nod_m                                    = y_nod_array(sa2_macropore_point_a.idx_nod);
sa2_macropore_point_a.z_nod_m                                    = z_nod_array(sa2_macropore_point_a.idx_nod);
sa2_macropore_point_a.legend                                     = sprintf('x = %i ; y = %i; z = %i',sa2_macropore_point_a.x,sa2_macropore_point_a.z);
sa2_macropore_point_a.mark_shape                                 = 'cx';

sa3_macropore_point_a.x                                          = 100;
[M index_sa3_macropore_point_a]                                  = min(abs(sa3_macropore_point_a.x-x_array) );
sa3_macropore_point_a.y                                          = -2 * pi * sa3_macropore_point_a.x;
sa3_macropore_point_a.z                                          = depth_groundwater_table_m + 1;
xdist_array                                                      = sa3_macropore_point_a.x-x_nod_array; ydist_array = sa3_macropore_point_a.y-y_nod_array;zdist_array = sa3_macropore_point_a.z-z_nod_array;
dist_array                                                       = (xdist_array .^ 2+ydist_array .^ 2+zdist_array .^ 2).^0.5;
[~,sa3_macropore_point_a.idx_nod]                                = min(dist_array);
sa3_macropore_point_a.x_nod_m                                    = x_nod_array(sa3_macropore_point_a.idx_nod);
sa3_macropore_point_a.y_nod_m                                    = y_nod_array(sa3_macropore_point_a.idx_nod);
sa3_macropore_point_a.z_nod_m                                    = z_nod_array(sa3_macropore_point_a.idx_nod);
sa3_macropore_point_a.legend                                     = sprintf('x = %i ; y = %i; z = %i',sa3_macropore_point_a.x,sa3_macropore_point_a.z);
sa3_macropore_point_a.mark_shape                                 = 'cx';

sa1_macropore_point_a.x                                          = 140;
[M index_sa1_macropore_point_a]                                  = min(abs(sa1_macropore_point_a.x-x_array) );
sa1_macropore_point_a.y                                          = -2 * pi * sa1_macropore_point_a.x;
sa1_macropore_point_a.z                                          = depth_groundwater_table_m + 1;
xdist_array                                                      = sa1_macropore_point_a.x-x_nod_array; ydist_array = sa1_macropore_point_a.y-y_nod_array;zdist_array = sa1_macropore_point_a.z-z_nod_array;
dist_array                                                       = (xdist_array .^ 2+ydist_array .^ 2+zdist_array .^ 2).^0.5;
[~,sa1_macropore_point_a.idx_nod]                                = min(dist_array);
sa1_macropore_point_a.x_nod_m                                    = x_nod_array(sa1_macropore_point_a.idx_nod);
sa1_macropore_point_a.y_nod_m                                    = y_nod_array(sa1_macropore_point_a.idx_nod);
sa1_macropore_point_a.z_nod_m                                    = z_nod_array(sa1_macropore_point_a.idx_nod);
sa1_macropore_point_a.legend                                     = sprintf('x = %i ; y = %i; z = %i',sa1_macropore_point_a.x,sa1_macropore_point_a.z);
sa1_macropore_point_a.mark_shape                                 = 'cx';


%%
vertical_line_pond_centre.x                                      = 1.0;
abs_xdist_array_m                                                = abs(vertical_line_pond_centre.x-x_nod_array);
min_abs_xdist                                                    = min(abs_xdist_array_m);
vertical_line_pond_centre.idx_nod                                = find( abs_xdist_array_m == min_abs_xdist);
vertical_line_pond_centre.idx_nod_pore                           = vertical_line_pond_centre.idx_nod(1:inp.nn1);
vertical_line_pond_centre.idx_nod_middle                         = vertical_line_pond_centre.idx_nod(inp.nn1 + 1:2*inp.nn1);
vertical_line_pond_centre.idx_nod_matrix                         = vertical_line_pond_centre.idx_nod(2*inp.nn1 + 1:end);
abs_xdist_array_m                                                = abs(vertical_line_pond_centre.x-x_ele_array);
min_abs_xdist                                                    = min(abs_xdist_array_m);
vertical_line_pond_centre.idx_ele                                = find( abs_xdist_array_m == min_abs_xdist);
vertical_line_pond_centre.idx_ele_pore                           = vertical_line_pond_centre.idx_ele(1:inp.nn1-1);
vertical_line_pond_centre.idx_ele_matrix                         = vertical_line_pond_centre.idx_ele(inp.nn1:end);

vertical_line_pond_boundary.x                                    = pond_radius_m-20;    %50
abs_xdist_array_m                                                = abs(vertical_line_pond_boundary.x-x_nod_array);
min_abs_xdist                                                    = min(abs_xdist_array_m);
vertical_line_pond_boundary.idx_nod                              = find( abs_xdist_array_m == min_abs_xdist);
vertical_line_pond_boundary.idx_nod_pore                         = vertical_line_pond_boundary.idx_nod(1:inp.nn1);
vertical_line_pond_boundary.idx_nod_middle                       = vertical_line_pond_boundary.idx_nod(inp.nn1 + 1:2*inp.nn1);
vertical_line_pond_boundary.idx_nod_matrix                       = vertical_line_pond_boundary.idx_nod(2*inp.nn1 + 1:end);
abs_xdist_array_m                                                = abs(vertical_line_pond_boundary.x-x_ele_array);
min_abs_xdist                                                    = min(abs_xdist_array_m);
vertical_line_pond_boundary.idx_ele                              = find( abs_xdist_array_m == min_abs_xdist);
vertical_line_pond_boundary.idx_ele_pore                         = vertical_line_pond_boundary.idx_ele(1:inp.nn1-1);
vertical_line_pond_boundary.idx_ele_matrix                       = vertical_line_pond_boundary.idx_ele(inp.nn1:end);

vertical_line_pond_outside.x                                     = pond_radius_m+20;    %50
abs_xdist_array_m                                                = abs(vertical_line_pond_outside.x-x_nod_array);
min_abs_xdist                                                    = min(abs_xdist_array_m);
vertical_line_pond_outside.idx_nod                               = find( abs_xdist_array_m == min_abs_xdist);
vertical_line_pond_outside.idx_nod_pore                          = vertical_line_pond_outside.idx_nod(1:inp.nn1);
vertical_line_pond_outside.idx_nod_middle                        = vertical_line_pond_outside.idx_nod(inp.nn1 + 1:2*inp.nn1);
vertical_line_pond_outside.idx_nod_matrix                        = vertical_line_pond_outside.idx_nod(2*inp.nn1 + 1:end);
abs_xdist_array_m                                                = abs(vertical_line_pond_outside.x-x_ele_array);
min_abs_xdist                                                    = min(abs_xdist_array_m);
vertical_line_pond_outside.idx_ele                               = find( abs_xdist_array_m == min_abs_xdist);
vertical_line_pond_outside.idx_ele_pore                          = vertical_line_pond_outside.idx_ele(1:inp.nn1-1);
vertical_line_pond_outside.idx_ele_matrix                        = vertical_line_pond_outside.idx_ele(inp.nn1:end);


fig_pos.left                                                     = 0.05;
fig_pos.bottom                                                   = 0.83;
fig_pos.length                                                   = 0.13;
fig_pos.height                                                   = 0.08;
fig_pos.xy_plot_length                                           = 0.08;
fig_pos.length_legend                                            = 0.05;
%nt                                                              =10;
%a.fig                                                           =figure;
a.fs                                                             = 8;
a.fs_legend                                                      = 5;
a.lw                                                             = 1.5;
a.cz                                                             = 8;
fs                                                               = 5; % sampling frequency
% Creating the movie : quality                                   = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt                                                               = 90;
%set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]);  % maximize the plotting figure

%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig                                                            = figure;
set(gcf,'Units','normalized', 'OuterPosition',[0 0 1.0 1.0]);  % maximize the plotting figure
mov                                                              =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
mov.FrameRate                                                    = 2;mov.Quality=qt;
open(mov);

for nt = 1:200:length(ele)

     %% -------- contour plot on saturation in the macropore ---------
    a.sub11  = subplot('position'...
         ,[fig_pos.left,fig_pos.bottom,...
          fig_pos.length,fig_pos.height]);
    
    % write pressure and conc in matrix form.
    s_matrix                                                     = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2,inp.nn3]);
    contourf(x_array,z_array,squeeze(s_matrix(:,1,:) ) );hold on;
    plot(x_nod_mtx_3d(source_point_index(1:ny:end) ),z_nod_mtx_3d(source_point_index(1:ny:end) ),'b*','MarkerSize',a.cz-2);
    plot(x_nod_mtx_3d(sink_point_index(1:ny:end) ),z_nod_mtx_3d(sink_point_index(1:ny:end) ),'r*','MarkerSize',a.cz-2);
    xline(sa2_macropore_point_a.x,'r-','linewidth',a.lw);
    xline(sa3_macropore_point_a.x,'g-','linewidth',a.lw);   
    xline(sa1_macropore_point_a.x,'b-','linewidth',a.lw);       
    colorbar('eastoutside');
    get(gca,'xtick');
    set(gca,'linewidth',a.lw); 
    grid on;    
    set(gca,'fontsize',a.fs,'fontweight','bold');
    title('Saturation in macropores','FontSize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
	caxis([0 1]);hold off;
    %axis([10, 40,9,10])
    %% -------- contour plot on concentration in the macropore ---------
    a.sub21    = subplot('position'...
         ,[fig_pos.left+fig_pos.length+0.23*fig_pos.length,fig_pos.bottom,...
          fig_pos.length,fig_pos.height]);    
    % write pressure and conc in matrix form.
    c_matrix                                                     = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2,inp.nn3]);
    contourf(x_array,z_array,squeeze(c_matrix(:,1,:) ) );
    colorbar('eastoutside');
    get(gca,'xtick');
    set(gca,'linewidth',a.lw); 
    grid on;    
    set(gca,'fontsize',a.fs,'fontweight','bold');
    title('Concentration in macropores','FontSize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);    
    ylabel('Elevation (m)','FontSize',a.fs);
	caxis([0 0.01]);
    %axis([10, 40,9,10])    
%      %% -------- plot velocity in the macropore  ---------    
%     a.sub31                                                    = subplot('position'...
%          ,[fig_pos.left+2*fig_pos.length+0.5*fig_pos.length,fig_pos.bottom,...
%           fig_pos.length,fig_pos.height]);    
%     vx_array                                                   = cell2mat(arrayfun(@(y) y.terms{vx_idx},ele(nt),'UniformOutput',false) );
%     vy_array                                                   = cell2mat(arrayfun(@(y) y.terms{vy_idx},ele(nt),'UniformOutput',false) );
%     vz_array                                                   = cell2mat(arrayfun(@(y) y.terms{vz_idx},ele(nt),'UniformOutput',false) );
%     vx_mtx                                                     = reshape(vx_array,[nez, ney, nex]);
%     vy_mtx                                                     = reshape(vy_array,[nez, ney, nex]);
%     vz_mtx                                                     = reshape(vz_array,[nez, ney, nex]);    
%     a.plot31                                                   = contourf(squeeze(x_ele_mtx_3d(:,1,:) ),squeeze(z_ele_mtx_3d(:,1,:) ),squeeze(vy_mtx(:,1,:) ) .* squeeze( porosity_ele_mtx(:,1,:) ) * c.secPday);hold on;
% 	set(gca,'ColorScale','log');
%     colorbar;caxis([-1e-20 1]);
%     a.plot31                                                   = quiver(squeeze(x_ele_mtx_3d(:,1,:) ),squeeze(z_ele_mtx_3d(:,1,:) ),squeeze(vx_mtx(:,1,:) ) .* squeeze( porosity_ele_mtx(:,1,:) ) * c.secPday,squeeze(vz_mtx(:,1,:) ) .* squeeze( porosity_ele_mtx(:,1,:) ) * c.secPday, 'color','r' );
%     set(a.plot31,'AutoScale','on', 'AutoScaleFactor',0.5)
%     hold off
%     grid on;
%     xlim([0 max(x_array)]);
%     ylim([min(z_array) max(z_array)]);    
%     set(gca,'fontsize',a.fs,'fontweight','bold');
%     title('Velocity (m/day) in macropore','FontSize',a.fs,'fontweight','bold');
%     xlabel('Distance (m)','FontSize',a.fs);    
%     ylabel('Elevation (m)','FontSize',a.fs);   
    %% -------------  sub 4 Evaporation and saturation over time  --------------    
    a.sub51    = subplot('position'...
         ,[fig_pos.left+2*fig_pos.length+0.45*fig_pos.length,...
           fig_pos.bottom,...
           fig_pos.length,fig_pos.height]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
	[Ax,Line1,Line2]                                             = plotyy( x_nod_array(sink_point_index(1:ny:end) ), cell2mat(arrayfun(@(y) y.terms{s_idx}(sink_point_index(1:ny:end) ),nod(nt),'UniformOutput',false) ),...
                              x_nod_array(sink_point_index(1:ny:end) ), -et_mmPday(1:ny:end,nt) );	
%     hold(Ax(2) );
% 	ax                                                           = plot(Ax(2),time_measured_data_day(1:100:end),...
%                   pet_movmean_mmPday(1:100:end),...
%                   'r.');hold off 

%     a.plot3                                                    = plotyy([time_measured_data_day(1:1000:end),time_nod_day(1:nt)],...
%                    [Tbl.pet_mmPday(1:1000:end),-et_mmPday(1:nt)],....
%                   time_nod_day(1:nt)',...
%                   Saturation(1:nt,sink_point_index_gravity_compensated(1) )'...
%                     );         
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    set(Ax,'FontSize',a.fs,'fontweight','bold');
	ylim([0 10]);
    set(Ax(1),'XLim',[0 max(x_array(:) )],'YLim',[0 1.1],'YTick',[0, 0.5, 1],'fontsize',a.fs,'linewidth',a.lw);
	set(Ax(2),'XLim',[0 max(x_array(:) )],'YLim',[0 10],'YTick',[0, 2, 4, 6],'fontsize',a.fs,'linewidth',a.lw);
	set(Line1,'linewidth',a.lw);
	set(Line2,'linewidth',a.lw);	
    set(gca,'linewidth',a.lw); 
    grid on;
%	set(Ax(3),'XLim',[0 180],'YLim',[0 10]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);
	ylabel(Ax(1),{'Saturation'},'FontSize',a.fs);    % right y-axis   
    ylabel(Ax(2),{'ET (mm/day)'},'FontSize',a.fs,'fontweight','bold'); % left y-axis
    hold off;
                
%% --------- sub 5-1 plot saturation profiles in macropore
    a.sub51            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+0.8*fig_pos.length,...
         fig_pos.bottom,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);
    plot(nod(nt).terms{s_idx}(vertical_line_pond_centre.idx_nod_pore)...
         ,z_array...
         ,'r-','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{s_idx}(vertical_line_pond_boundary.idx_nod_pore)...
         ,z_array...
         ,'g-','linewidth',a.lw) ;
    plot(nod(nt).terms{s_idx}(vertical_line_pond_outside.idx_nod_pore)...
         ,z_array...
         ,'b-','linewidth',a.lw) ;
        
    set(gca,'linewidth',a.lw); 
    grid on;
    %plot(slab(2,:,1),slab(1,:,1),'rd',slab(2,:,2),...
    %    slab(1,:,2),'go',slab(2,:,3),slab(1,:,3),'cx','linewidth',a.lw);hold off;
    %axis([-0.7*fig_pos.height, 1.05,0,0.05])
    x_rectangle                                                  = [-0.2,1.2,1.2,-0.2];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');hold off;
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');hold off
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
    xlim([-0.2 1.2]);
% 	set(p3,'FaceAlpha',0.5);
    txt  = sprintf('Result at day %.2f , k_silt                     = %.2e m2, k_macropore= %.2e m2, k_sand = %.2e m2, \n k_pond= %.2e m2,por_silt= %.2e,por_macropore= %.2e,por_pond= %.2e, init P in clay layer = %.2e Pa'...
         ,nod(nt).tout*c.dayPsec,permeability_silt_m2,permeability_macropore_m2,permeability_sand_m2,permeability_pond_m2,porosity_silt,porosity_macropore,porosity_pond,initial_silt_layer_pressure_pa);  
    title(txt,'FontSize',a.fs-2);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Saturation (-)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    
    %% --------- sub 5-2 plot conc profiles  ------------------
    a.sub52            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+fig_pos.xy_plot_length+0.95*fig_pos.length,...
         fig_pos.bottom,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);

    plot(nod(nt).terms{c_idx}(vertical_line_pond_centre.idx_nod_pore)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_centre.idx_nod_pore)...
         ,'r-','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{c_idx}(vertical_line_pond_boundary.idx_nod_pore)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_boundary.idx_nod_pore)...
         ,'g-','linewidth',a.lw) ;
    plot(nod(nt).terms{c_idx}(vertical_line_pond_outside.idx_nod_pore)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_outside.idx_nod_pore)...
         ,'b-','linewidth',a.lw) ;    
    x_rectangle                                                  = [0,0.04,0.04,0];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');
    hold off;
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
% 	set(p3,'FaceAlpha',0.5);
    set(gca,'linewidth',a.lw); 
    grid on;
    %rectangle('Position',[0 -5 1 5],'FaceColor','red','FaceAlpha',0.5);  % no facealpha
    %plot(clab(2,:,1),clab(1,:,1),'rd',clab(2,:,2),clab(1,:,2)...
    %,'go',clab(2,:,3),clab(1,:,3),'cx','MarkerSize',a.cz,'linewidth',a.lw);hold off
    %axis([-0.7*fig_pos.height, 0.28,0,0.05])
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Salinity (kg/kg)','FontSize',a.fs);
%     ylabel('Elevation (m)','FontSize',a.fs);
    
    %% --------- sub 5-3 plot x and y velocity  ------------------
    a.sub53            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+2*fig_pos.xy_plot_length + 1.1*fig_pos.length,...
         fig_pos.bottom,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);
    plot(ele(nt).terms{vx_idx}(1:inp.nn1-1) * c.secPday .* inp.por(1:inp.nn1-1)...
         ,ele(nt).terms{zele_idx}(1:inp.nn1-1)...
         ,'linewidth',a.lw) ;hold on
    plot(ele(nt).terms{vz_idx}(1:inp.nn1-1) * c.secPday .* inp.por(1:inp.nn1-1) ...
         ,ele(nt).terms{zele_idx}(1:inp.nn1-1)...
         ,'linewidth',a.lw); 
    hold off
    x_rectangle                                                  = [-0.5,0.5,0.5,-0.5];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
% 	set(p3,'FaceAlpha',0.5);
    set(gca,'linewidth',a.lw); 
    grid on;    
    %rectangle('Position',[0 -5 1 5],'FaceColor','red','FaceAlpha',0.5);  % no facealpha
    %plot(clab(2,:,1),clab(1,:,1),'rd',clab(2,:,2),clab(1,:,2)...
    %,'go',clab(2,:,3),clab(1,:,3),'cx','MarkerSize',a.cz,'linewidth',a.lw);hold off
    %axis([-0.7*fig_pos.height, 0.28,0,0.05])
    get(gca,'xtick');
    xlim([-0.5 0.5]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Vx(blue) Vz(red) (m/day)','FontSize',a.fs);
%     ylabel('Elevation (m)','FontSize',a.fs);   
    %% --------- sub 5-4 plot relative permeability  ------------------
    a.sub54            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+3*fig_pos.xy_plot_length + 1.4*fig_pos.length,...
         fig_pos.bottom,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);
    semilogx(ele(nt).terms{rek_idx}(vertical_line_pond_centre.idx_ele_pore) ...
         ,ele(nt).terms{zele_idx}(vertical_line_pond_centre.idx_ele_pore)...
         ,'r-','linewidth',a.lw) ;hold on
    semilogx(ele(nt).terms{rek_idx}(vertical_line_pond_boundary.idx_ele_pore) ...
         ,ele(nt).terms{zele_idx}(vertical_line_pond_boundary.idx_ele_pore)...
         ,'g-','linewidth',a.lw) ; 
    semilogx(ele(nt).terms{rek_idx}(vertical_line_pond_outside.idx_ele_pore) ...
         ,ele(nt).terms{zele_idx}(vertical_line_pond_outside.idx_ele_pore)...
         ,'b-','linewidth',a.lw) ;     
    x_rectangle                                                  = [1e-20,1,1,1e-20];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
% 	set(p3,'FaceAlpha',0.5);
    set(gca,'linewidth',a.lw); 
    grid on;    
    get(gca,'xtick');
    xlim([1e-20 1]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Relative permeability (-)','FontSize',a.fs);hold off;
%     ylabel('Elevation (m)','FontSize',a.fs);    
    
          %% -------- contour plot on saturation in the middle layer ---------
    a.sub12    = subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-fig_pos.height-0.75*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    % write pressure and conc in matrix form.
    s_matrix                                                     = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2,inp.nn3]);
    contourf(x_array,z_array,squeeze(s_matrix(:,2,:) ) );hold on;
    plot(x_nod_mtx_3d(source_point_index(2:ny:end) ),z_nod_mtx_3d(source_point_index(2:ny:end) ),'b*','MarkerSize',a.cz-2);
    plot(x_nod_mtx_3d(sink_point_index(2:ny:end) ),z_nod_mtx_3d(sink_point_index(2:ny:end) ),'r*','MarkerSize',a.cz-2);
    colorbar('eastoutside');
    get(gca,'xtick');
    set(gca,'linewidth',a.lw); 
    grid on;    
    set(gca,'fontsize',a.fs,'fontweight','bold');
    title('Saturation in the middle layer','FontSize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
	caxis([0 1]);hold off;
    %axis([10, 40,9,10])

     %% -------- contour plot on concentration in the middle layer ---------
    a.sub22            = subplot('position'...
         ,[fig_pos.left+fig_pos.length+0.23*fig_pos.length,fig_pos.bottom-fig_pos.height-0.75*fig_pos.height,...
          fig_pos.length,fig_pos.height]);   
    % write pressure and conc in matrix form.
    c_matrix                                                     = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2,inp.nn3]);
    contourf(x_array,z_array,squeeze(c_matrix(:,2,:) ) );
    colorbar('eastoutside');
    get(gca,'xtick');
    set(gca,'linewidth',a.lw); 
    grid on;    
    set(gca,'fontsize',a.fs,'fontweight','bold');
    title('Concentration in the middle layer','FontSize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
	caxis([0 0.01]);
    %axis([10, 40,9,10])
  
    %% -------------  sub 4 Evaporation and saturation in the middle layer over time  --------------    
    a.sub51            = subplot('position'...
         ,[fig_pos.left+2*fig_pos.length+0.45*fig_pos.length,...
           fig_pos.bottom-fig_pos.height-0.75*fig_pos.height,...
           fig_pos.length,fig_pos.height]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
	[Ax,Line1,Line2]                                             = plotyy( x_nod_array(sink_point_index(2:ny:end) ), cell2mat(arrayfun(@(y) y.terms{s_idx}(sink_point_index(2:ny:end) ),nod(nt),'UniformOutput',false) ),...
                              x_nod_array(sink_point_index(2:ny:end) ), -et_mmPday(2:ny:end,nt) );
    hold(Ax(2) );
% 	ax                                                           = plot(Ax(2),time_measured_data_day(1:100:end),...
%                   pet_movmean_mmPday(1:100:end),...
%                   'r.');hold off 
%     a.plot3                                                    = plot(Ax(2),time_nod_day(1:nt),arrayfun(@(y) y.terms{s_idx}(sink_point_index(6) ),nod(1:nt) ),'g-');
%     a.plot3                                                    = plot(Ax(2),time_nod_day(1:nt),arrayfun(@(y) y.terms{s_idx}(sink_point_index(9) ),nod(1:nt) ),'b-' );

%     a.plot3                                                    = plotyy([time_measured_data_day(1:1000:end),time_nod_day(1:nt)],...
%                    [Tbl.pet_mmPday(1:1000:end),-et_mmPday(1:nt)],....
%                   time_nod_day(1:nt)',...
%                   Saturation(1:nt,sink_point_index_gravity_compensated(1) )'...
%                     );         
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    set(Ax,'FontSize',a.fs,'fontweight','bold');
	ylim([0 10]);
    set(Ax(1),'XLim',[0 max(x_array(:) )],'YLim',[0 1.1],'YTick',[0, 0.5, 1],'fontsize',a.fs,'linewidth',a.lw);
	set(Ax(2),'XLim',[0 max(x_array(:) )],'YLim',[0 10],'YTick',[0, 2, 4, 6],'fontsize',a.fs,'linewidth',a.lw);
	set(Line1,'linewidth',a.lw);
	set(Line2,'linewidth',a.lw);	
    set(gca,'linewidth',a.lw); 
    grid on;
%	set(Ax(3),'XLim',[0 180],'YLim',[0 10]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);
	ylabel(Ax(1),{'Saturation'},'FontSize',a.fs);    % right y-axis   
    ylabel(Ax(2),{'ET (mm/day)'},'FontSize',a.fs,'fontweight','bold'); % left y-axis   ho   
    hold off;                                
       %% --------- sub 6-1 plot saturation profiles in the middle layer
    a.sub61            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+0.8*fig_pos.length,...
         fig_pos.bottom-fig_pos.height-0.75*fig_pos.height,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);
    
    plot(nod(nt).terms{s_idx}(vertical_line_pond_centre.idx_nod_pore+nez)...
         ,z_array...
         ,'r--','linewidth',a.lw) ;
    plot(nod(nt).terms{s_idx}(vertical_line_pond_boundary.idx_nod_pore+nez)...
         ,z_array...
         ,'g--','linewidth',a.lw) ;
    plot(nod(nt).terms{s_idx}(vertical_line_pond_outside.idx_nod_pore+nez)...
         ,z_array...
         ,'b--','linewidth',a.lw) ;   
     
    set(gca,'linewidth',a.lw); 
    grid on;
    %plot(slab(2,:,1),slab(1,:,1),'rd',slab(2,:,2),...
    %    slab(1,:,2),'go',slab(2,:,3),slab(1,:,3),'cx','linewidth',a.lw);hold off;
    %axis([-0.7*fig_pos.height, 1.05,0,0.05])
    x_rectangle                                                  = [-0.2,1.2,1.2,-0.2];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');hold off;
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');hold off
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
    xlim([-0.2 1.2]);
% 	set(p3,'FaceAlpha',0.5);	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Saturation (-)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    %% --------- sub 6-2 plot conc profiles in the middle layer  ------------------
    a.sub62            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+fig_pos.xy_plot_length+0.95*fig_pos.length,...
         fig_pos.bottom-fig_pos.height-0.75*fig_pos.height,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);

    plot(nod(nt).terms{c_idx}(vertical_line_pond_centre.idx_nod_pore+nz)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_centre.idx_nod_pore+nz)...
         ,'r--','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{c_idx}(vertical_line_pond_boundary.idx_nod_pore+nz)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_boundary.idx_nod_pore+nz)...
         ,'g--','linewidth',a.lw) ;
    plot(nod(nt).terms{c_idx}(vertical_line_pond_outside.idx_nod_pore+nz)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_outside.idx_nod_pore+nz)...
         ,'b--','linewidth',a.lw) ;    
    x_rectangle                                                  = [0,0.04,0.04,0];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');
    hold off;
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
% 	set(p3,'FaceAlpha',0.5);
    set(gca,'linewidth',a.lw); 
    grid on;
    %rectangle('Position',[0 -5 1 5],'FaceColor','red','FaceAlpha',0.5);  % no facealpha
    %plot(clab(2,:,1),clab(1,:,1),'rd',clab(2,:,2),clab(1,:,2)...
    %,'go',clab(2,:,3),clab(1,:,3),'cx','MarkerSize',a.cz,'linewidth',a.lw);hold off
    %axis([-0.7*fig_pos.height, 0.28,0,0.05])
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Salinity (kg/kg)','FontSize',a.fs);
%     ylabel('Elevation (m)','FontSize',a.fs);

      %% -------- contour plot on saturation in the matrix ---------
    a.sub13            = subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    % write pressure and conc in matrix form.
    s_matrix                                                     = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2,inp.nn3]);
    contourf(x_array,z_array,squeeze(s_matrix(:,3,:) ) );hold on;
    plot(x_nod_mtx_3d(source_point_index(3:ny:end) ),z_nod_mtx_3d(source_point_index(3:ny:end) ),'b*','MarkerSize',a.cz-2);
    plot(x_nod_mtx_3d(sink_point_index(3:ny:end) ),z_nod_mtx_3d(sink_point_index(3:ny:end) ),'r*','MarkerSize',a.cz-2);
    colorbar('eastoutside');
    get(gca,'xtick');
    set(gca,'linewidth',a.lw); 
    grid on;    
    set(gca,'fontsize',a.fs,'fontweight','bold');
    title('Saturation in the matrix','FontSize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);    
    ylabel('Elevation (m)','FontSize',a.fs);
	caxis([0 1]);hold off;
    %axis([10, 40,9,10])


     %% -------- contour plot on concentration in the matrix  ---------
    a.sub22            = subplot('position'...
         ,[fig_pos.left+fig_pos.length+0.23*fig_pos.length,fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
          fig_pos.length,fig_pos.height]);    
    % write pressure and conc in matrix form.
    c_matrix                                                     = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2,inp.nn3]);
    contourf(x_array,z_array,squeeze(c_matrix(:,3,:) ) );
    colorbar('eastoutside');
    get(gca,'xtick');
    set(gca,'linewidth',a.lw); 
    grid on;    
    set(gca,'linewidth',a.lw); 
    set(gca,'fontsize',a.fs,'fontweight','bold');
    title('Concentration in the matrix','FontSize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
	caxis([0 0.01]);
    %axis([10, 40,9,10])  
    

%          %% -------- plot velocity in the Matrix  ---------    
%     a.sub32                                                    = subplot('position'...
%          ,[fig_pos.left+2*fig_pos.length+0.5*fig_pos.length,fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
%           fig_pos.length,fig_pos.height]);    
%     vx_array                                                   = cell2mat(arrayfun(@(y) y.terms{vx_idx},ele(nt),'UniformOutput',false) );
%     vy_array                                                   = cell2mat(arrayfun(@(y) y.terms{vy_idx},ele(nt),'UniformOutput',false) );
%     vz_array                                                   = cell2mat(arrayfun(@(y) y.terms{vz_idx},ele(nt),'UniformOutput',false) );
%     vx_mtx                                                     = reshape(vx_array,[nez, ney, nex]);
%     vy_mtx                                                     = reshape(vy_array,[nez, ney, nex]);
%     vz_mtx                                                     = reshape(vz_array,[nez, ney, nex]);    
%     a.plot32                                                   = contourf(squeeze(x_ele_mtx_3d(:,2,:) ),squeeze(z_ele_mtx_3d(:,2,:) ),squeeze(vy_mtx(:,2,:) ) .*  squeeze( porosity_ele_mtx(:,2,:) ) * c.ms2mmday);hold on;
% 	set(gca,'ColorScale','log');   
% 	colorbar;caxis([-1e-20 11]);
%     a.plot32                                                   = quiver(squeeze(x_ele_mtx_3d(:,2,:) ),squeeze(z_ele_mtx_3d(:,2,:) ),squeeze(vx_mtx(:,2,:) ) .*  squeeze( porosity_ele_mtx(:,2,:) ) *c.ms2mmday,squeeze(vz_mtx(:,2,:) ) .*  squeeze( porosity_ele_mtx(:,2,:) ) *c.ms2mmday, 'color','r');
%     set(a.plot32,'AutoScale','on', 'AutoScaleFactor',0.5)
%     hold off
%     grid on; 
%     xlim([0 max(x_array)]);
%     ylim([min(z_array) max(z_array)]);      
%     set(gca,'fontsize',a.fs,'fontweight','bold');
%     title('Velocity (m/day) in the matrix','FontSize',a.fs,'fontweight','bold');
%     xlabel('Distance (m)','FontSize',a.fs);    
%     ylabel('Elevation (m)','FontSize',a.fs);
    %% -------------  sub 4 Evaporation and saturation in the matrix over time  --------------    
    a.sub51            = subplot('position'...
         ,[fig_pos.left+2*fig_pos.length+0.45*fig_pos.length,...
           fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
           fig_pos.length,fig_pos.height]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
	[Ax,Line1,Line2]                                             = plotyy( x_nod_array(sink_point_index(3:ny:end) ), cell2mat(arrayfun(@(y) y.terms{s_idx}(sink_point_index(3:ny:end) ),nod(nt),'UniformOutput',false) ),...
                              x_nod_array(sink_point_index(3:ny:end) ), -et_mmPday(3:ny:end,nt) );

%     a.plot3                                                    = plotyy([time_measured_data_day(1:1000:end),time_nod_day(1:nt)],...
%                    [Tbl.pet_mmPday(1:1000:end),-et_mmPday(1:nt)],....
%                   time_nod_day(1:nt)',...
%                   Saturation(1:nt,sink_point_index_gravity_compensated(1) )'...
%                     );         
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    set(Ax,'FontSize',a.fs,'fontweight','bold');
	ylim([0 10]);
    set(Ax(1),'XLim',[0 max(x_array(:) )],'YLim',[0 1.1],'YTick',[0, 0.5, 1],'fontsize',a.fs,'linewidth',a.lw);
	set(Ax(2),'XLim',[0 max(x_array(:) )],'YLim',[0 10],'YTick',[0, 2, 4, 6],'fontsize',a.fs,'linewidth',a.lw);
	set(Line1,'linewidth',a.lw);
	set(Line2,'linewidth',a.lw);	
    set(gca,'linewidth',a.lw); 
    grid on;
%	set(Ax(3),'XLim',[0 180],'YLim',[0 10]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Distance (m)','FontSize',a.fs);
	ylabel(Ax(1),{'Saturation'},'FontSize',a.fs);    % right y-axis   
    ylabel(Ax(2),{'ET (mm/day)'},'FontSize',a.fs,'fontweight','bold'); % left y-axis   ho   
    hold off;            
       %% --------- sub 7-1 plot saturation profiles in matrix
    a.sub71            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+0.8*fig_pos.length,...
         fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);
    plot(nod(nt).terms{s_idx}(vertical_line_pond_centre.idx_nod_matrix)...
         ,z_array...
         ,'r:','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{s_idx}(vertical_line_pond_boundary.idx_nod_matrix)...
         ,z_array...
         ,'g:','linewidth',a.lw) ;
    plot(nod(nt).terms{s_idx}(vertical_line_pond_outside.idx_nod_matrix)...
         ,z_array...
         ,'b:','linewidth',a.lw) ;
        
    set(gca,'linewidth',a.lw); 
    grid on;
    %plot(slab(2,:,1),slab(1,:,1),'rd',slab(2,:,2),...
    %    slab(1,:,2),'go',slab(2,:,3),slab(1,:,3),'cx','linewidth',a.lw);hold off;
    %axis([-0.7*fig_pos.height, 1.05,0,0.05])
    x_rectangle                                                  = [-0.2,1.2,1.2,-0.2];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');hold off;
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');hold off
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
    xlim([-0.2 1.2]);
% 	set(p3,'FaceAlpha',0.5);	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Saturation (-)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    
    %% --------- sub 7-2 plot conc profiles  ------------------
    a.sub72            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+fig_pos.xy_plot_length+0.95*fig_pos.length,...
         fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);

    plot(nod(nt).terms{c_idx}(vertical_line_pond_centre.idx_nod_matrix)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_centre.idx_nod_matrix)...
         ,'r:','linewidth',a.lw) ;hold on
    plot(nod(nt).terms{c_idx}(vertical_line_pond_boundary.idx_nod_matrix)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_boundary.idx_nod_matrix)...
         ,'g:','linewidth',a.lw) ;
    plot(nod(nt).terms{c_idx}(vertical_line_pond_outside.idx_nod_matrix)...
         ,nod(nt).terms{z_idx}(vertical_line_pond_outside.idx_nod_matrix)...
         ,'b:','linewidth',a.lw) ;    
    x_rectangle                                                  = [0,0.04,0.04,0];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');
    hold off;
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
% 	set(p3,'FaceAlpha',0.5);
    set(gca,'linewidth',a.lw); 
    grid on;
    %rectangle('Position',[0 -5 1 5],'FaceColor','red','FaceAlpha',0.5);  % no facealpha
    %plot(clab(2,:,1),clab(1,:,1),'rd',clab(2,:,2),clab(1,:,2)...
    %,'go',clab(2,:,3),clab(1,:,3),'cx','MarkerSize',a.cz,'linewidth',a.lw);hold off
    %axis([-0.7*fig_pos.height, 0.28,0,0.05])
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Salinity (kg/kg)','FontSize',a.fs);
%     ylabel('Elevation (m)','FontSize',a.fs);
    
    %% --------- sub 7-3 plot x and y velocity  ------------------
    a.sub73            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+2*fig_pos.xy_plot_length + 1.1*fig_pos.length,...
         fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);
    plot(ele(nt).terms{vx_idx}(inp.nn1:2*(inp.nn1-1) ) * c.secPday .* inp.por(inp.nn1:2*(inp.nn1-1) )...
         ,ele(nt).terms{zele_idx}(1:inp.nn1-1)...
         ,'r:','linewidth',a.lw) ;hold on
    plot(ele(nt).terms{vz_idx}(inp.nn1:2*(inp.nn1-1) ) * c.secPday .* inp.por(inp.nn1:2*(inp.nn1-1) ) ...
         ,ele(nt).terms{zele_idx}(1:inp.nn1-1)...
         ,'b:','linewidth',a.lw); 
    hold off
    x_rectangle                                                  = [-0.5,0.5,0.5,-0.5];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
% 	set(p3,'FaceAlpha',0.3);
    set(gca,'linewidth',a.lw); 
    grid on;    
    %rectangle('Position',[0 -5 1 5],'FaceColor','red','FaceAlpha',0.5);  % no facealpha
    %plot(clab(2,:,1),clab(1,:,1),'rd',clab(2,:,2),clab(1,:,2)...
    %,'go',clab(2,:,3),clab(1,:,3),'cx','MarkerSize',a.cz,'linewidth',a.lw);hold off
    %axis([-0.7*fig_pos.height, 0.28,0,0.05])
    get(gca,'xtick');
    xlim([-0.5 0.5]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Vx(blue) Vz(red) (m/day)','FontSize',a.fs);
%     ylabel('Elevation (m)','FontSize',a.fs);   

    %% --------- sub 7-4 plot relative permeability velocity  ------------------
    a.sub74            = subplot('position'...
         ,[fig_pos.left+3*fig_pos.length+3*fig_pos.xy_plot_length + 1.4*fig_pos.length,...
         fig_pos.bottom-2*fig_pos.height-1.5*fig_pos.height,...
          fig_pos.xy_plot_length,fig_pos.height+0.01]);
    semilogx(ele(nt).terms{rek_idx}(vertical_line_pond_centre.idx_ele_matrix) ...
         ,ele(nt).terms{zele_idx}(vertical_line_pond_centre.idx_ele_matrix)...
         ,'r:','linewidth',a.lw) ;hold on
    semilogx(ele(nt).terms{rek_idx}(vertical_line_pond_boundary.idx_ele_matrix) ...
         ,ele(nt).terms{zele_idx}(vertical_line_pond_boundary.idx_ele_matrix)...
         ,'g:','linewidth',a.lw) ; 
    semilogx(ele(nt).terms{rek_idx}(vertical_line_pond_outside.idx_ele_matrix) ...
         ,ele(nt).terms{zele_idx}(vertical_line_pond_outside.idx_ele_matrix)...
         ,'b:','linewidth',a.lw) ;     
    x_rectangle                                                  = [1e-20,1,1,1e-20];
    z_rectangle1                                                 = [depth_groundwater_table_m,depth_groundwater_table_m,-10,-10];
    z_rectangle2                                                 = [1,1,-1,-1];
% 	z_rectangle3                                                 = [-1,-1,-2,-2];
    p1                                                           = patch(x_rectangle,z_rectangle1,'k');
	p2                                                           = patch(x_rectangle,z_rectangle2,'k');
% 	p3                                                           = patch(x_rectangle,z_rectangle3,'b');
    set(p1,'FaceAlpha',0.3);
	set(p2,'FaceAlpha',0.3);
% 	set(p3,'FaceAlpha',0.5);
    set(gca,'linewidth',a.lw); 
    grid on;    
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    xlabel('Relative permeability (-)','FontSize',a.fs); 
    xlim([1e-20 1.2]);hold off;

        %% -------------   Groundwater table over time outside the basin SA1 --------------     
    a.sub41            = subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-3*fig_pos.height-2.2*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    a.plot41                                                     = plot(time_measured_data_day(1:1000:end),...
																	Tbl.sa1_groundwater_table_rise_mm(1:1000:end)*c.mPmm,'k.','displayname','Mea. at SA1');    
    hold on; 
    a.plot41                                                     = plot(time_measured_data_day(1:1000:end),...
																	Tbl.sa4_groundwater_table_rise_mm(1:1000:end)*c.mPmm,'k:','displayname','Mea. at SA4');        
    %a.plot11                                                    = plot(time_nod_day(1:nt),watertable_diff_sa1_m(1:nt),...
    %         'r-','linewidth',a.lw);
%     a.plot41                                                   = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod),nod(1:nt) ) - nod(1).terms{p_idx}(sa1_macropore_point_a.idx_nod + 1) )  /c.rhow_pure_water/c.g ,'r-','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k-','linewidth',a.lw,'displayname','simulation'   );    
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r-','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod + 1),nod(2) ) )  /c.   rhow_pure_water/c.g ,'g-','linewidth',a.lw,'displayname','simulation'   );
	a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2),nod(2) ) )  /c.				rhow_pure_water/c.g ,'b-','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+nz-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k--','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+nz),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r--','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+nz + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+nz + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g--','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+nz+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+nz+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b--','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2*nz-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k:','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2*nz),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r:','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2*nz + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2*nz + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g:','linewidth',a.lw,'displayname','simulation'   );
    a.plot41                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2*nz+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa1_macropore_point_a.idx_nod+2*nz+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b:','linewidth',a.lw,'displayname','simulation'   );
     
    set(gca,'linewidth',a.lw); 
    grid on;
    ylim([-3 5.5]);
    yticks([-2 1 2.5]);
    xlim([0 floor(time_nod_day(end) )]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
%     legend('fontsize',a.fs_legend,'Location',[fig_pos.length+0.2*fig_pos.length,fig_pos.bottom,fig_pos.length_legend,fig_pos.height]);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel({'Groundwater';'Level';'at SA1 (m)'},'FontSize',a.fs);   hold off

            %% -------------   Groundwater table over time within the basin SA2 --------------      
    a.sub42            = subplot('position'...
         ,[fig_pos.left+fig_pos.length+0.45*fig_pos.length, fig_pos.bottom-3*fig_pos.height-2.2*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    a.plot42                                                     = plot(time_measured_data_day(1:1000:end),...
            Tbl.sa2_groundwater_table_rise_mm(1:1000:end)*c.mPmm,'k.','displayname','Mea. at SA1');hold on  
    %a.plot11                                                    = plot(time_nod_day(1:nt),watertable_diff_sa1_m(1:nt),...
    %         'r-','linewidth',a.lw);
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k-','linewidth',a.lw,'displayname','simulation'   );    
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r-','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g-','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b-','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+nz-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k--','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+nz),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r--','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+nz + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+nz + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g--','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+nz+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+nz+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b--','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2*nz-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k:','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2*nz),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r:','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2*nz + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2*nz + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g:','linewidth',a.lw,'displayname','simulation'   );
    a.plot42                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2*nz+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa2_macropore_point_a.idx_nod+2*nz+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b:','linewidth',a.lw,'displayname','simulation'   );

    set(gca,'linewidth',a.lw); 
    grid on;
    ylim([-3 6.0]);
    yticks([-2 1 2.5]);
    xlim([0 floor(time_nod_day(end) )]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
%     legend('fontsize',a.fs_legend,'Location',[fig_pos.length+0.2*fig_pos.length,fig_pos.bottom,fig_pos.length_legend,fig_pos.height]);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel({'Groundwater';'Level';'at SA2 (m)'},'FontSize',a.fs); hold off    
            %% -------------   Groundwater table over time within the basin SA3 --------------      
    a.sub43            = subplot('position'...
         ,[fig_pos.left+2*fig_pos.length+0.9*fig_pos.length, fig_pos.bottom-3*fig_pos.height-2.2*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    a.plot43                                                     = plot(time_measured_data_day(1:1000:end),...
            Tbl.sa3_groundwater_table_rise_mm(1:1000:end)*c.mPmm,'k.','displayname','Mea. at SA1');hold on  
    %a.plot11                                                    = plot(time_nod_day(1:nt),watertable_diff_sa1_m(1:nt),...
    %         'r-','linewidth',a.lw);
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k-','linewidth',a.lw,'displayname','simulation'   );    
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r-','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g-','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b-','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+nz-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k--','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+nz),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r--','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+nz + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+nz + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g--','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+nz+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+nz+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b--','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2*nz-1),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'k:','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2*nz),nod(1:nt) ) - 0)  /c.rhow_pure_water/c.g ,'r:','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2*nz + 1),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2*nz + 1),nod(2) ) )  /c.rhow_pure_water/c.g ,'g:','linewidth',a.lw,'displayname','simulation'   );
    a.plot43                                                     = plot(time_nod_day(1:nt) , (arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2*nz+2),nod(1:nt) ) - arrayfun(@(y) y.terms{p_idx}(sa3_macropore_point_a.idx_nod+2*nz+2),nod(2) ) )  /c.rhow_pure_water/c.g ,'b:','linewidth',a.lw,'displayname','simulation'   );
    set(gca,'linewidth',a.lw); 
    grid on;
    ylim([-3 6.0]);
    yticks([-2 1 2.5]);
    xlim([0 floor(time_nod_day(end) )]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
%     legend('fontsize',a.fs_legend,'Location',[fig_pos.length+0.2*fig_pos.length,fig_pos.bottom,fig_pos.length_legend,fig_pos.height]);
    xlabel('Time (day)','FontSize',a.fs);   
    ylabel({'Groundwater';'Level';'at SA3 (m)'},'FontSize',a.fs);  hold off   


%     vy_matrix                                                  = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1,inp.nn3-1]);
%     [depth_m depth_index]=min(z_ele_mtx_3d(:,1,1)-(-5) );
%     contourf(squeeze(x_ele_mtx_3d(depth_index,:,:) ),squeeze(y_ele_mtx_3d(depth_index,:,:) ),squeeze(vy_matrix(depth_index,:,:) ) );
%     colorbar('eastoutside');
%     caxis([-1e-6 1e-6]);
%     get(gca,'xtick');
%     set(gca,'linewidth',a.lw); 
%     grid on;    
%     set(gca,'fontsize',a.fs,'fontweight','bold');
%     title('Velocity in Y direction (m/s)','FontSize',a.fs,'fontweight','bold');
%     xlabel('X (m)','FontSize',a.fs);
%     ylabel('Y (m)','FontSize',a.fs);    
         %% -------------  surface water level over time  --------------
    a.sub311                                                     =subplot('position'...
         ,[fig_pos.left+3*fig_pos.length + 1.35*fig_pos.length , ...
           fig_pos.bottom-3*fig_pos.height-2.2*fig_pos.height , ...
           fig_pos.length , fig_pos.height]);         
%     [M,I]                                                      = min(abs(time_nod_day(nt)-time_measured_data_day) )  ;
    a.plot311                                                    = plot(time_measured_data_day(1:1000:end),...
        (Tbl.p3_cs451(1:1000:end)-0.2213),'b.');hold on   
    pond_depth_cell_integrated_m_sa2_macropore                   =   squeeze(sum(water_height_per_cell_mtx_m(1:index_location_of_pond_boundary_nod(index_point_1_pond),1,index_point_1_pond,:) ) );
    pond_depth_cell_integrated_m_sa2_mid                         =   squeeze(sum(water_height_per_cell_mtx_m(1:index_location_of_pond_boundary_nod(index_point_1_pond),2,index_point_1_pond,:) ) );
    pond_depth_cell_integrated_m_sa2_matrix                      =   squeeze(sum(water_height_per_cell_mtx_m(1:index_location_of_pond_boundary_nod(index_point_1_pond),3,index_point_1_pond,:) ) );
    


    %     pond_depth_cell_integrated_below_1_m_1                 =   squeeze(sum(water_height_per_cell_mtx_m(1:index_location_of_pond_boundary_nod(index_point_1_pond) + 1,3,index_point_1_pond,:) ) );
%     pond_depth_cell_integrated_below_2_m_1                     =   squeeze(sum(water_height_per_cell_mtx_m(1:index_location_of_pond_boundary_nod(index_point_1_pond)+2,3,index_point_1_pond,:) ) );
    a.plot311                                                    = plot(time_nod_day(1:nt) , pond_depth_cell_integrated_m_sa2_macropore(1:nt),'r-','linewidth',a.lw  );
    a.plot311                                                    = plot(time_nod_day(1:nt) , pond_depth_cell_integrated_m_sa2_mid(1:nt),'r--','linewidth',a.lw  );
    a.plot311                                                    = plot(time_nod_day(1:nt) , pond_depth_cell_integrated_m_sa2_matrix(1:nt),'r:','linewidth',a.lw  );

%     legend('Measurement','Macropore','Middle layer','Matrix','Position',[0.45 0.58 0.02 0.01]);
    %     a.plot32                                               = plot(time_nod_day(1:nt) , pond_depth_cell_integrated_below_1_m_1(1:nt),'g--','linewidth',a.lw  );
%     a.plot32                                                   = plot(time_nod_day(1:nt) , pond_depth_cell_integrated_below_2_m_1(1:nt),'g.','linewidth',a.lw  ); 
    hold off
    set(gca,'linewidth',a.lw); 
    grid on;
    xlim([0 floor(time_nod_day(end) )]);
    ylim([0 1.2]);
	yticks([0.6 1.2]);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel({'Surface water';'depth';'(m)'},'FontSize',a.fs); 
    set(gca,'fontsize',a.fs,'fontweight','bold');
        %% -------------  sub 4 Evaporation and saturation at one sink point in macropore over time  --------------    
    a.sub51            = subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-4*fig_pos.height-2.6*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
    [Ax,Line1,Line2]                                             = plotyy(time_nod_day(1:nt),[(arrayfun(@(y) y.terms{s_idx}(sink_point_index(1) ),nod(1:nt) ) )',...
                                                  (arrayfun(@(y) y.terms{s_idx}(sink_point_index(4) ),nod(1:nt) ) )',...
                                                  (arrayfun(@(y) y.terms{s_idx}(sink_point_index(7) ),nod(1:nt) ) )'],...
													time_nod_day(1:nt),[-et_mmPday(1,1:nt)',...
													-et_mmPday(4,1:nt)',...
													-et_mmPday(7,1:nt)']); 
    set(Line1(1),'color','b','linestyle','-','linewidth',a.lw);
    set(Line1(2),'color','b','linestyle','--','linewidth',a.lw);
    set(Line1(3),'color','b','linestyle',':','linewidth',a.lw);

    set(Line2(1),'color','r','linestyle','-','linewidth',a.lw);
    set(Line2(2),'color','r','linestyle','--','linewidth',a.lw);
    set(Line2(3),'color','r','linestyle',':','linewidth',a.lw);

%     hold(Ax);

%     a.plot3                                                    = plotyy([time_measured_data_day(1:1000:end),time_nod_day(1:nt)],...
%                    [Tbl.pet_mmPday(1:1000:end),-et_mmPday(1:nt)],....
%                   time_nod_day(1:nt)',...
%                   Saturation(1:nt,sink_point_index_gravity_compensated(1) )'...
%                     );         
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    set(Ax,'FontSize',a.fs,'fontweight','bold');
	ylim([0 10]);
    set(Ax(1),'XLim',[0 floor(time_nod_day(end) )],'YLim',[0 1.1],'YTick',[0, 0.5, 1],'fontsize',a.fs,'linewidth',a.lw);
	set(Ax(2),'XLim',[0 floor(time_nod_day(end) )],'YLim',[0 10],'YTick',[0, 2, 4, 6],'fontsize',a.fs,'linewidth',a.lw);
	set(Line1,'linewidth',a.lw);
	set(Line2,'linewidth',a.lw);	
    set(gca,'linewidth',a.lw); 
    grid on;
%	set(Ax(3),'XLim',[0 180],'YLim',[0 10]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
%     title('Saturation and ET in macropores',a.fs,'fontweight','bold');
    xlabel('Time (day)','FontSize',a.fs);
	ylabel(Ax(1),{'Saturation'},'FontSize',a.fs);    % right y-axis   
    ylabel(Ax(2),{'ET - macropore';'(mm/day)'},'FontSize',a.fs,'fontweight','bold');hold off; % left y-axis   
       %% -------------  sub 4 Evaporation and saturation at one sink point in middle over time  --------------    
    a.sub51            = subplot('position'...
         ,[fig_pos.left+fig_pos.length+0.45*fig_pos.length ,fig_pos.bottom-4*fig_pos.height-2.6*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
    [Ax,Line1,Line2]                                             = plotyy(time_nod_day(1:nt),[(arrayfun(@(y) y.terms{s_idx}(sink_point_index(2) ),nod(1:nt) ) )',...
                                                  (arrayfun(@(y) y.terms{s_idx}(sink_point_index(5) ),nod(1:nt) ) )',...
                                                  (arrayfun(@(y) y.terms{s_idx}(sink_point_index(8) ),nod(1:nt) ) )'],...
													time_nod_day(1:nt),[-et_mmPday(2,1:nt)',...
													-et_mmPday(5,1:nt)',...
													-et_mmPday(8,1:nt)']);
    set(Line1(1),'color','b','linestyle','-','linewidth',a.lw);
    set(Line1(2),'color','b','linestyle','--','linewidth',a.lw);
    set(Line1(3),'color','b','linestyle',':','linewidth',a.lw);

    set(Line2(1),'color','r','linestyle','-','linewidth',a.lw);
    set(Line2(2),'color','r','linestyle','--','linewidth',a.lw);
    set(Line2(3),'color','r','linestyle',':','linewidth',a.lw);
                                              
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    set(Ax,'FontSize',a.fs,'fontweight','bold');
	ylim([0 10]);
    set(Ax(1),'XLim',[0 floor(time_nod_day(end) )],'YLim',[0 1.1],'YTick',[0, 0.5, 1],'fontsize',a.fs,'linewidth',a.lw);
	set(Ax(2),'XLim',[0 floor(time_nod_day(end) )],'YLim',[0 10],'YTick',[0, 2, 4, 6],'fontsize',a.fs,'linewidth',a.lw);
	set(Line1,'linewidth',a.lw);
	set(Line2,'linewidth',a.lw);	
    set(gca,'linewidth',a.lw); 
    grid on;
%	set(Ax(3),'XLim',[0 180],'YLim',[0 10]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
%     title('Saturation and ET in the middle',a.fs,'fontweight','bold');
    xlabel('Time (day)','FontSize',a.fs);
	ylabel(Ax(1),{'Saturation'},'FontSize',a.fs);    % right y-axis   
    ylabel(Ax(2),{'ET - middle';'(mm/day)'},'FontSize',a.fs,'fontweight','bold');hold off; % left y-axis      
       %% -------------  sub 4 Evaporation and saturation at one sink point in matrix over time  --------------    
    a.sub51            = subplot('position'...
         ,[fig_pos.left+2*fig_pos.length+0.9*fig_pos.length ,fig_pos.bottom-4*fig_pos.height-2.6*fig_pos.height,...
          fig_pos.length,fig_pos.height]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
    %arrayfun(@(y) y.tout,nod(1:nt) )
	%ylim(ax,[0 10]);
    [Ax,Line1,Line2]                                             = plotyy(time_nod_day(1:nt),[(arrayfun(@(y) y.terms{s_idx}(sink_point_index(3) ),nod(1:nt) ) )',...
                                                  (arrayfun(@(y) y.terms{s_idx}(sink_point_index(6) ),nod(1:nt) ) )',...
                                                  (arrayfun(@(y) y.terms{s_idx}(sink_point_index(9) ),nod(1:nt) ) )'],...
													time_nod_day(1:nt),[-et_mmPday(3,1:nt)',...
													-et_mmPday(6,1:nt)',...
													-et_mmPday(9,1:nt)']);
%     a.plot3                                                    = plotyy([time_measured_data_day(1:1000:end),time_nod_day(1:nt)],...
%                    [Tbl.pet_mmPday(1:1000:end),-et_mmPday(1:nt)],....
%                   time_nod_day(1:nt)',...
%                   Saturation(1:nt,sink_point_index_gravity_compensated(1) )'...
%                     );         
    set(Line1(1),'color','b','linestyle','-','linewidth',a.lw);
    set(Line1(2),'color','b','linestyle','--','linewidth',a.lw);
    set(Line1(3),'color','b','linestyle',':','linewidth',a.lw);

    set(Line2(1),'color','r','linestyle','-','linewidth',a.lw);
    set(Line2(2),'color','r','linestyle','--','linewidth',a.lw);
    set(Line2(3),'color','r','linestyle',':','linewidth',a.lw);

    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'fontweight','bold');
    set(Ax,'FontSize',a.fs,'fontweight','bold');
	ylim([0 10]);
    set(Ax(1),'XLim',[0 floor(time_nod_day(end) )],'YLim',[0 1.1],'YTick',[0, 0.5, 1],'fontsize',a.fs,'linewidth',a.lw);
	set(Ax(2),'XLim',[0 floor(time_nod_day(end) )],'YLim',[0 10],'YTick',[0, 2, 4, 6],'fontsize',a.fs,'linewidth',a.lw);
	set(Line1,'linewidth',a.lw);
	set(Line2,'linewidth',a.lw);	
    set(gca,'linewidth',a.lw); 
    grid on;
%	set(Ax(3),'XLim',[0 180],'YLim',[0 10]);
%     set(gca,'XTicklabel',[]);
    set(gca,'fontsize',a.fs,'fontweight','bold');
%     title('Saturation and ET in the matrix',a.fs,'fontweight','bold');    
    xlabel('Time (day)','FontSize',a.fs);
	ylabel(Ax(1),{'Saturation'},'FontSize',a.fs);    % right y-axis   
    ylabel(Ax(2),{'ET - matrix';'(mm/day)'},'FontSize',a.fs,'fontweight','bold');hold off; % left y-axis          
        %% -------------  Mass balance over time  --------------   
    a.sub32            = subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-5.3*fig_pos.height-4.3*fig_pos.height...
          1.2*fig_pos.length,2.5*fig_pos.height]);         
    saturation_ele_mtx                                           = nod_to_ele(nex,ney,nez,s_matrix);
    vz_matrix                                                    = reshape(ele(nt).terms{vz_idx},[inp.nn1-1,inp.nn2-1,inp.nn3-1]);
    a.plot32                                                     = plot(time_nod_day(1:nt) , input_total_kg(1:nt)/c.mlPl,'r-','linewidth',a.lw,'DisplayName','Input' );hold on;
    h11                                                          = plot(time_measured_data_day,Tbl.pump+Tbl.rainfall_volume_ML,'r--','linewidth',a.lw,'DisplayName','Input-obs');
    %     ylim([0 120]);
    %     plot(time_nod_day(1:nt),water_outflow_right_boundary_kg(1:nt)/c.mlPl,'g-','linewidth',a.lw);
    h2                                                           = plot(time_nod_day(1:nt),outflow_bcop_cumsum_kg(1:nt)/c.mlPl,'b-','linewidth',a.lw,'DisplayName','Outflow from boundary');
    %     plot(time_nod_day(1:nt),water_increased_kg(1:nt)/1e6,'k-','linewidth',a.lw);
    %     plot(time_nod_day(1:nt),(water_outflow_right_boundary_kg(1:nt)+water_increased_kg(1:nt) )/1e6,'y-','linewidth',a.lw);
    h3                                                           = plot(time_nod_day(1:nt),(pond_water_sum_kg(1:nt)-pond_water_sum_kg(1) )/c.mlPl,'m-','linewidth',a.lw,'DisplayName','Water of pond');
    h33                                                          = plot(time_measured_data_day,Tbl.volumeTOTALML,'m--','linewidth',a.lw,'DisplayName','Water of pond-obs');
%     h4                                                         = plot(time_nod_day(1:nt),water_outflow_from_pond_boundary_xz_kg(1:nt)/c.mlPl,'c-','linewidth',a.lw);
%     h5                                                         = plot(time_nod_day(1:nt),water_flux_groundwater_table_xz_ele_cumsum(1:nt)/c.mlPl,'linewidth',a.lw);
    h6                                                           = plot(time_nod_day(1:nt),cumsum(sum(-c.rhow_pure_water*area_xy_nod_mtx_m2_3d(sink_point_index) .* et_mmPday(:,1:nt)*c.mPmm*c.dayPsec .* time_series(1:nt),1) )/c.mlPl,...
                                        'color','#0A952C','LineStyle','-','linewidth',a.lw,'DisplayName','ET');
    h66                                                          = plot(time_measured_data_day,pet_observed_kg/c.mlPl,'color','#0A952C','LineStyle','--','linewidth',a.lw,'DisplayName','ET-obs');
    h7                                                           = plot(time_nod_day(1:nt),(Infiltration_kg(1:nt)-Infiltration_kg(1) )/c.mlPl,'k-','linewidth',a.lw,'DisplayName','Infiltration');%%infiltration from model
    h77                                                          = plot(time_measured_data_day,Tbl.infiltration_cum_ML,'k--','linewidth',a.lw,'DisplayName','Infiltration-obs');    
    h8                                                           = plot(time_nod_day(1:nt),silt_zone_water_cell_kg(1:nt)/c.mlPl,'linewidth',a.lw,'DisplayName','Water in the clay layer');
%     h9                                                         = plot(time_nod_day(1:nt),input_total_kg(1:nt)/c.mlPl-pond_water_sum_kg(1:nt)/c.mlPl-...
%                cumsum(sum(-c.rhow_pure_water*area_xy_mtx_m2(sink_point_index_gravity_compensated) .* et_mmPday(:,1:nt)*c.mPmm*c.dayPsec .* time_series(1:nt),1) )/c.mlPl-silt_zone_water_cell_xyt_kg(1:nt)/c.mlPl,'linewidth',a.lw);%%infiltration from model
    %     a.plot1                                                = plot(time_measured_data_day(1:1000:end), Tbl.sa2_groundwater_table_rise_mm(1:1000:end)*c.mPmm,'b.','displayname','Mea. at SA2' );hold on
    hold off

    ylabel({'Water';'Volume';'ML'},'FontSize',a.fs);
    xlabel('Time (day)')
    xlim([0 floor(time_nod_day(end) )]);
    set(gca,'fontsize',a.fs+2,'fontweight','bold');
    set(gca,'linewidth',a.lw); 
%     set(get(get(h11,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(h33,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(h66,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     set(get(get(h77,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	legend('position',[fig_pos.left + 1.4*fig_pos.length,...
                       fig_pos.bottom-5*fig_pos.height-3.6*fig_pos.height,...
                       fig_pos.length,fig_pos.height]);
    grid on;
    set(gca,'linewidth',a.lw); 
    grid on;
    xlim([0 floor(time_nod_day(end) )]);
    ylim([0 130]);
	yticks([50 120]);
    ylabel({'Mass (ML)'},'FontSize',a.fs); 
    xlabel('Time (day)','FontSize',a.fs);
    set(gca,'fontsize',a.fs,'fontweight','bold');        
        %% -------------  Mass exchange cross layers over time  --------------   
    a.sub32            = subplot('position'...
         ,[fig_pos.left+2.9*fig_pos.length,fig_pos.bottom-5.3*fig_pos.height-4.3*fig_pos.height...
          1.2*fig_pos.length,2.5*fig_pos.height]);  
%     h11                                                        = plot ( time_nod_day( 1:nt ), nod_to )
    h1                                                           = plot(time_nod_day(1:nt),exchange_macropore_kg(1:nt) / c.mlPl,'linewidth',a.lw,'DisplayName','Exchange - macropores');hold on;
    h2                                                           = plot(time_nod_day(1:nt),exchange_matrix_kg(1:nt) / c.mlPl,'linewidth',a.lw,'DisplayName','Exchange - matrix');
    h3                                                           = plot(time_nod_day(1:nt),(Infiltration_macropore_kg(1:nt) - Infiltration_macropore_kg(1) ) / c.mlPl,'linewidth',a.lw,'DisplayName','Infiltration - macropores');
    h4                                                           = plot(time_nod_day(1:nt),(Infiltration_matrix_kg(1:nt) - Infiltration_matrix_kg(1) ) / c.mlPl,'linewidth',a.lw,'DisplayName','Infiltration - matrix');
    h5                                                           = plot(time_nod_day(1:nt),(Recharge_macropore_kg(1:nt) - Recharge_macropore_kg(1) ) / c.mlPl,'linewidth',a.lw,'DisplayName','Recharge - macropores');
    h6                                                           = plot(time_nod_day(1:nt),(Recharge_matrix_kg(1:nt) - Recharge_matrix_kg(1) ) / c.mlPl,'linewidth',a.lw,'DisplayName','Recharge - matrix');

    legend('position',[fig_pos.left+4.4*fig_pos.length,...
                       fig_pos.bottom-4*fig_pos.height-4.5*fig_pos.height,...
                       0.5*fig_pos.length,fig_pos.height]);
    grid on;
    set(gca,'linewidth',a.lw); 
    grid on;
    xlim([0 floor(time_nod_day(end) )]);
    ylim([0 max(Infiltration_macropore_kg)/c.mlPl]);
	yticks([0 0.5*max((Infiltration_macropore_kg - Infiltration_macropore_kg(1) ) )/c.mlPl max((Infiltration_macropore_kg - Infiltration_macropore_kg(1) ) )/c.mlPl]);
    ylabel({'Mass (ML)'},'FontSize',a.fs); 
    xlabel('Time (day)','FontSize',a.fs);  
    set(gca,'fontsize',a.fs,'fontweight','bold');        
    hold off;
         %% -------------  flux at the bottom  --------------   
   
    a.sub8             = subplot('position'...
                    ,[fig_pos.left+5.7*fig_pos.length,fig_pos.bottom-3*fig_pos.height-2.2*fig_pos.height...
                    fig_pos.length,fig_pos.height]);        
    h1                                                           = plot(squeeze(x_ele_mtx_3d(1,1,:) ),squeeze(sum(aquitard_flux_ele_xy_kgPs(:,1,:,nt),1) )/c.dayPsec,'r-','linewidth',a.lw);hold on;
    h2                                                           = plot(squeeze(x_ele_mtx_3d(1,1,:) ),squeeze(sum(aquitard_flux_ele_xy_kgPs(:,2,:,nt),1) )/c.dayPsec,'r:','linewidth',a.lw);
    xlim([0 max(squeeze(x_ele_mtx_3d(1,1,:) ) )]);
%     ylim([min(aquitard_flux_ele_xy_kgPs(:) )/c.dayPsec max(aquitard_flux_ele_xy_kgPs(:) )/c.dayPsec]);
    xlabel('Distance (m)','FontSize',a.fs);      
    ylabel({'Flux (kg/day)'},'FontSize',a.fs); 
    set(gca,'fontsize',a.fs,'fontweight','bold');        
    set(gca,'linewidth',a.lw); 
    grid on;
    hold off;
         %% -------------  Distribution of BCOP points  --------------   
    a.sub91            = subplot('position'...
                ,[fig_pos.left+5.7*fig_pos.length,fig_pos.bottom-4*fig_pos.height-2.9*fig_pos.height...
                fig_pos.length,fig_pos.height]);    
    for n                                                        = 1:length(ipbc_node_idx_array)/ny        
        plot(time_nod_day(1:nt),cumsum(outflow_bcop_node_kg(n,1:nt) )/c.mlPl,'linewidth',a.lw);hold on;
    end
    hold off
    xlim([0 floor(time_nod_day(end) )]);
    ylim([0 2*cumsum(sum(outflow_bcop_node_kg(1,:),2) )/c.mlPl]);
    ylabel({'Mass from ','boundary (ML)'},'FontSize',a.fs); 
    xlabel('Time (day)','FontSize',a.fs);    
    set(gca,'fontsize',a.fs,'fontweight','bold');        
    set(gca,'linewidth',a.lw); 
    grid on;
   
    a.sub92            = subplot('position'...
                ,[fig_pos.left+5.7*fig_pos.length,fig_pos.bottom-5*fig_pos.height-3.5*fig_pos.height...
                fig_pos.length,fig_pos.height]);    
    for n                                                        = length(ipbc_node_idx_array)/ny + 1:2*length(ipbc_node_idx_array)/ny        
        plot(time_nod_day(1:nt),cumsum(outflow_bcop_node_kg(n,1:nt) )/c.mlPl,'linewidth',a.lw);hold on;
    end
    hold off
    xlim([0 floor(time_nod_day(end) )]);
    ylim([0 4*cumsum(sum(outflow_bcop_node_kg(1,:),2) )/c.mlPl]);
    ylabel({'Mass from','boundary (ML)'},'FontSize',a.fs); 
    xlabel('Time (day)','FontSize',a.fs); 
    set(gca,'fontsize',a.fs,'fontweight','bold');        
    set(gca,'linewidth',a.lw); 
    grid on;
    a.sub93            = subplot('position'...
                ,[fig_pos.left+5.7*fig_pos.length,fig_pos.bottom-6*fig_pos.height-4.0*fig_pos.height...
                fig_pos.length,fig_pos.height]);    
    for n                                                        = 2*length(ipbc_node_idx_array)/ny + 1:3*length(ipbc_node_idx_array)/ny        
        plot(time_nod_day(1:nt),cumsum(outflow_bcop_node_kg(n,1:nt) )/c.mlPl,'linewidth',a.lw);hold on;
    end
    hold off
    xlim([0 floor(time_nod_day(end) )]);
    ylim([0 2*cumsum(sum(outflow_bcop_node_kg(1,:),2) )/c.mlPl]);
    ylabel({'Mass from','boundary (ML)'},'FontSize',a.fs); 
    xlabel('Time (day)','FontSize',a.fs);        
    set(gca,'fontsize',a.fs,'fontweight','bold');        
    set(gca,'linewidth',a.lw); 
    grid on;

%          %% -------------  SWCC over time  --------------   
%    
%     a.sub8                                                     =subplot('position'...
%                     ,[fig_pos.left+5.5*fig_pos.length,fig_pos.bottom-5.3*fig_pos.height-4.3*fig_pos.height...
%                     0.5*fig_pos.length,2.5*fig_pos.height]);                
%     for n                                                      =1:swcc.no_session
%         hh                                                     =semilogy(swcc.por(n)*swcc.sw_fayer1995wrr{n},-swcc.psim,swcc.line_type{n},'displayname',swcc.soil_type{n});hold on
%     end 
%     xlabel('Vol. water content(-)','fontweight','bold','fontsize',a.fs)
%     ylabel('-Matric potential (m)','fontweight','bold','fontsize',a.fs)
%     grid on;
% %     legend ('fontsize',a.fs-1);
%     set(gca,'FontSize',a.fs,'FontWeight','bold','linewidth',2);hold off;
%          %% -------------  SWCC vs k over time  --------------   
%     
%     a.sub9                                                     =subplot('position'...
%                     ,[fig_pos.left+6.4*fig_pos.length,fig_pos.bottom-5.3*fig_pos.height-4.3*fig_pos.height...
%                     0.5*fig_pos.length,2.5*fig_pos.height]);                
%     for n                                                      =1:swcc.no_session
%         hh                                                     =semilogy(swcc.por(n)*swcc.sw_fayer1995wrr{n},swcc.ksat(n)*swcc.kr_mualem1976wrr{n},swcc.line_type{n},'displayname',swcc.soil_type{n});hold on
%     end
%     xlabel('Vol. water content(-)','fontweight','bold','fontsize',a.fs)
%     ylabel('Permeability (m^2)','fontweight','bold','fontsize',a.fs)
%     grid on;
% 	legend('position',[fig_pos.left+5.8*fig_pos.length,...
%                        fig_pos.bottom-3*fig_pos.height-3.6*fig_pos.height,...
%                        fig_pos.length,fig_pos.height]);hold off;
    set(gca,'FontSize',a.fs,'FontWeight','bold','linewidth',2);        
    saveas(gcf,['chart',num2str(nt),'.png']);
    F                                                            = getframe(gcf); % save the current figure
    writeVideo(mov,F);% add it as the next frame of the movie
end

close(mov);
close all;
close(mov);
close all;

