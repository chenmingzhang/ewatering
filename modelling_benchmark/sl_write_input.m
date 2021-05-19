clear all
close all
c=ConstantObj();  % all the constants, we suggest to use c.g rather than 9.81 in the script to enhance readerbility;

%% SUTRA.fil

fil=filObj('read_from_file','no');

% the default value for 'fid' is -1, meaning, this will not write to SUTRA.fil if export
fil.terms.('inp').('fname')  = 'PART1.inp' ; fil.terms.('inp').('fid')   = 50;
fil.terms.('ics').('fname')  = 'PART1.ics' ; fil.terms.('ics').('fid')   = 55;
fil.terms.('lst').('fname')  = 'PART1.lst' ; fil.terms.('lst').('fid')   = 60;
fil.terms.('rst').('fname')  = 'PART1.rst' ; fil.terms.('rst').('fid')   = 66;
fil.terms.('nod').('fname')  = 'PART1.nod' ; fil.terms.('nod').('fid')   = 30;
fil.terms.('ele').('fname')  = 'PART1.ele' ; fil.terms.('ele').('fid')   = 40;
fil.terms.('obs').('fname')  = 'PART1.obs' ; fil.terms.('obs').('fid')   = 70;
fil.terms.('smy').('fname')  = 'PART1.smy' ; fil.terms.('smy').('fid')   = 80;
fil.terms.('bcof').('fname') = 'PART1.bcof'; fil.terms.('bcof').('fid')  = 91;
fil.terms.('bcos').('fname') = 'PART1.bcos'; fil.terms.('bcos').('fid')  = 93;
fil.terms.('bcop').('fname') = 'PART1.bcop'; fil.terms.('bcop').('fid')  = 92;
fil.terms.('bcou').('fname') = 'PART1.bcou'; fil.terms.('bcou').('fid')  = 94;


fil.export_to_file();

c_saltwater_kgPkg          = 0.035;
c_freshwater_kgPkg         = 0.001;
initial_head_aquifer_m     = -2.8;
initial_pond_water_depth_m = 0.5;
permeability_silt_m2       = 5.38e-15;
permeability_sand_m2       = 5.38e-11;
pond_radius_m              = 50.1;

%% inp file
dx      = 1.0;
dy      = 1;  % vertical direction
dz      = 1;
x_array = 0:dx:100;
y_array = -10:dy:1;
nx      = length(x_array);
ny      = length(y_array);

nex  = length(x_array)-1;
ney  = length(y_array)-1;

[x_nod_mtx,y_nod_mtx]=meshgrid(x_array,y_array);

keynodes            = zeros(size(x_nod_mtx,1),size(x_nod_mtx,2)+1);
keynodes(:,2:end-1) = (x_nod_mtx(:,1:end-1)+x_nod_mtx(:,2:end))/2;
keynodes(:,1)       = x_nod_mtx(:,1);
keynodes(:,end)     = x_nod_mtx(:,end);
dx_cell_mtx         = diff(keynodes,1,2); % check inbuilding function diff

keynodes            = zeros(size(y_nod_mtx,1)+1,size(y_nod_mtx,2));
keynodes(2:end-1,:) = (y_nod_mtx(1:end-1,:)+y_nod_mtx(2:end,:))/2;
keynodes(1,:)       = y_nod_mtx(1,:);
keynodes(end,:)     = y_nod_mtx(end,:);
dy_cell_mtx         = diff(keynodes,1,1); % check inbuilding function diff


nn       = nx*ny;
ne       = (nx-1)*(ny-1);
sequence = 'yxz';

%dataset 14
ii = (1:nn)';
x_nod_array  = x_nod_mtx(:);
y_nod_array  = y_nod_mtx(:);


%node_index_mtx = flip(reshape(ii,ny,nx));  %note node_index_mtx(52) = 52
node_index_mtx = reshape(ii,ny,nx);  %note node_index_mtx(52) = 52

% the reason of having a gravity compensated matrix is that the node reflects a xy domain with gravity facing down. 
% the original matrix has the same shape as the one with gravity compensation but it is bottom-up

node_index_mtx_gravity_compensated = flip(node_index_mtx);   % this makes a matrix that maps the node position in a xy plane.
y_nod_mtx_gravity_compensated_m      = flip(y_nod_mtx);
x_nod_mtx_gravity_compensated_m      = flip(x_nod_mtx);
dx_cell_mtx_gravity_compensated    = flip(dx_cell_mtx);

idx = 1;

iin1 = zeros(ne,1);
iin2 = zeros(ne,1);
iin3 = zeros(ne,1);
iin4 = zeros(ne,1);
x_ele_array=zeros(ne,1);
y_ele_array=zeros(ne,1);

for j  = 1:nex
    for i = 1:ney
        iin1(idx) = node_index_mtx(i,j);
        iin2(idx) = node_index_mtx(i,j+1);
        iin3(idx) = node_index_mtx(i+1,j+1);
        iin4(idx) = node_index_mtx(i+1,j);
        x_ele_array(idx)= mean([x_nod_array(iin1(idx)),x_nod_array(iin2(idx)),x_nod_array(iin3(idx)),x_nod_array(iin4(idx))] )   ;  % use mean method to get the middle of the cell
        y_ele_array(idx)= mean([y_nod_array(iin1(idx)),y_nod_array(iin2(idx)),y_nod_array(iin3(idx)),y_nod_array(iin4(idx))] )   ;  % use mean method to get the middle of the cell
        idx       = idx+1;

    end
end

x_ele_mtx_m = reshape(x_ele_array, ny-1,nx-1);
y_ele_mtx_m = reshape(y_ele_array, ny-1,nx-1);


x_ele_mtx_gravity_compensated_m = flip(x_ele_mtx_m);
y_ele_mtx_gravity_compensated_m = flip(y_ele_mtx_m);


%% dataset 17
%iqcp = -node_index_mtx_gravity_compensated(1,:)' ; % top cell,negative means evaporation is included
%qinc = dx_cell_mtx_gravity_compensated(1,:)' *dz ;
%uinc = zeros(size(iqcp)) +1e-4;   % evaporating water concentration




%ipbc = node_index_mtx_gravity_compensated(end,:)'; % only the last column. note that the result needs to be in a column
%pbc  = zeros(size(ipbc)) - 7.37325e+03;
%ubc  = zeros(size(ipbc)) + 3.0e-03;


% PART1 is the name
inp = inpObj('PART1','read_from_file','no');   % setup a empty inpObj

% dataset 1
inp.title1 = 'E-watering project';
inp.title2 = 'Generating input using sutralab';

% dataset 2a
%'SUTRA VERSION 2.2 SOLUTE TRANSPORT'
%'2D REGULAR MESH' 105 4001

inp.vermin = '2.2';  % note this should be a string not a number;
inp.simula = 'SOLUTE';


% dataset 2b
 %        2d mesh          ==>   ktype(1) = 2
 %        3d mesh          ==>   ktype(1) = 3
 %        irregular mesh   ==>   ktype(2) = 0
 %        layered mesh     ==>   ktype(2) = 1
 %        regular mesh     ==>   ktype(2) = 2
 %        blockwise mesh   ==>   ktype(2) = 3

inp.ktype(1)  = 2;  % 2D mesh
inp.mshtyp{1} = '2D';
inp.mshtyp{2} = 'REGULAR';

inp.nn1 = ny;
inp.nn2 = nx;


% ##  DATASET 3:  Simulation Control Numbers
inp.nn   = nn;
inp.ne   = ne;
inp.npbc = 0;   %length(ipbc);  revised after dataset 19
inp.nubc = 0;
inp.nsop = 0;       % dataset 17
inp.nsou = 0;
inp.nobs = 0;


%%##  DATASET 4:  Simulation Mode Options

inp.cunsat = 'UNSATURATED';
inp.cssflo = 'TRANSIENT FLOW';
inp.csstra = 'TRANSIENT TRANSPORT';
inp.cread  = 'COLD' ;
inp.istore = 9999;


%%##  DATASET 5:  Numerical Control Parameters
inp.up   = 0;
inp.gnup = 0.01;
inp.gnuu = 0.01;


%  DATASET 6:  Temporal Control and Solution Cycling Data
%  
inp.nsch  = 1;
inp.npcyc = 1;
inp.nucyc = 1;

%DATASET 6:  Temporal Control and Solution Cycling Data

inp.schnam = 'TIME_STEPS';
inp.schtyp = 'TIME CYCLE';
inp.creft  = 'ELAPSED';
%inp.scalt  = 6000;   %reduce to 3000
%inp.scalt  = 3000;   %reduce to 3000
inp.scalt  = 600;   %reduce to 3000
%inp.scalt  = 150;   %reduce to 3000
%inp.scalt  = 6000;   %reduce to 3000
%inp.scalt  = 0.6;   %reduce to 3000
inp.ntmax  = 10000;
inp.timei  = 0;
inp.timel  = 1.e99;
inp.timec  = 1.;
inp.ntcyc  = 9999;
inp.tcmult = 1;
inp.tcmin  = 1.e-20;
inp.tcmax  = 1.e99;


%##  DATASET 7:  ITERATION AND MATRIX SOLVER CONTROLS
%##  [ITRMAX]        [RPMAX]        [RUMAX]
inp.itrmax = 1000;
inp.rpmax  = 1e+5;
%inp.rpmax = 5e-2;  TO200317 too strigent for the first step
inp.rumax  = 1.0e-1;
% ##  [CSOLVP]  [ITRMXP]         [TOLP]
inp.csolvp = 'ORTHOMIN' ;
inp.itrmxp = 3000;
inp.tolp   = 1e-12;

%##  [CSOLVU]  [ITRMXU]         [TOLU]
inp.csolvu = 'ORTHOMIN';
inp.itrmxu = 3000;
inp.tolu   = 1e-12;


%##  DATASET 8:  Output Controls and Options
%## [NPRINT]  [CNODAL]  [CELMNT]  [CINCID]  [CPANDS]  [CVEL]  [CCORT] [CBUDG]   [CSCRN]  [CPAUSE]
%   2920        'N'        'N'        'N'        'Y'     'Y'        'Y'    'Y'      'Y' 'Y' 'Data Set 8A'

inp.nprint = 1000;
inp.cnodal = 'N';
inp.celmnt = 'N';
inp.cincid = 'N';
inp.cpands = 'Y';
inp.cvel   = 'Y';
inp.ccort  = 'N';
inp.cbudg  = 'Y';
inp.cscrn  = 'N';   % screen output
inp.cpause = 'Y';


%## [NCOLPR]    [NCOL]
%     -1000  'N'  'X'  'Y'  'P'  'U'  'S'  '-' 
%## [LCOLPR]    [LCOL]
%     1000 'E'  'X'  'Y'  'VX' 'VY' '-' 
%## [NOBCYC]    [INOB]
%
%##  [NBCFPR]  [NBCSPR]  [NBCPPR]  [NBCUPR]  [CINACT]
%     1000         1000       1000      1000       'Y'
%##
%##  DATASET 9:  Fluid Properties
%##     [COMPFL]           [CW]       [SIGMAW]        [RHOW0]       [URHOW0]        [DRWDU]        [VISC0]
%         0.0                1.         3.890D-10         1.0E+3             0.     7.0224E+02         1.0E-3
%##
%##  DATASET 10:  Solid Matrix Properties
%##     [COMPMA]           [CS]       [SIGMAS]         [RHOS]
%         0.0             0.             0.             1.
%##
%##  DATASET 11:  Adsorption Parameters


inp.ncolpr = -500;
%inp.ncol  = 'N'  'X'  'Y'  'P'  'U'  'S'  '-';
inp.ncol   = {['N'],['X' ],[ 'Y'  ],['P' ],[ 'U' ],[ 'S' ],[ '-']};

inp.lcolpr = 500;
inp.lcol   = {[ 'E' ],[ 'X' ],[ 'Y'  ],['VX' ],['VY' ],['-']};


inp.nbcfpr = 500;
inp.nbcspr = 500;
inp.nbcppr = 500;
inp.nbcupr = 500;
inp.cinact = 'Y';



%##    DATASET 9:  FLUID PROPERTIES
inp.compfl = 1.e-11;
inp.cw     = 1.;
inp.sigmaw = 1.e-9;
inp.rhow0  = 1000;
inp.urhow0 = 0;
inp.drwdu  = 700;
inp.visc0  = 0.001;


%##      DATASET 10:  SOLID MATRIX PROPERTIES
inp.compma = 1.e-7;
inp.cs     = 0;
inp.sigmas = 0;
inp.rhos   = 2600;   %solid density of sodium chloride

%##  DATASET 11:  Adsorption Parameters
%##     [ADSMOD]         [CHI1]         [CHI2]
%#'NONE'
%'FREUNDLICH' 1.D-47 0.05  #less rigid

%inp.adsmod = 'FREUNDLICH';
inp.adsmod = 'NONE';
inp.chi1   = 1.D-46;
inp.chi2   = 0.05 ;



%##
%##  DATASET 12:  Production of Energy or Solute Mass
%##     [PRODF0]       [PRODS0]       [PRODF1]       [PRODS1]
%0.             0.             0.             0.
inp.prodf0 = 0;
inp.prods0 = 0;
inp.prodf1 = 0;
inp.prods1 = 0;


%##
%##  DATASET 13:  Orientation of Coordinates to Gravity
%##      [GRAVX]        [GRAVY]        [GRAVZ]
%0.           -9.81          0.
inp.gravx = 0;
inp.gravy = -9.81;
inp.gravz = 0;

%##  DATASET 14:  NODEWISE DATA
%%##                              [SCALX] [SCALY] [SCALZ] [PORFAC]
inp.scalx  = 1.;
inp.scaly  = 1.;
inp.scalz  = 1.;
inp.porfac = 1.;

%##      [II]    [NRE    G(II)]  [X(II)] [Y(II)] [Z(II)] [POR(II)]
inp.ii   = (1:nn)';
inp.nreg = zeros(nn,1)+1;
inp.x    = x_nod_array;
inp.y    = y_nod_array;
%inp.z    = zeros(nn,1)+dz;
inp.z    = 2*pi*inp.x+0.1;
inp.por  = zeros(nn,1)+0.43;


%##                              [PMAXFA]        [PMINFA]        [ANG1FA]        [ALMAXF]        [ALMINF]        [ATMAXF]        [ATMINF]
%'ELEMENT'               1.0000000D+00   1.0000000D+00   1.0000000D+00   2 2 2 2
%##     [L]      [LREG(L)]       [PMAX(L)]       [PMIN(L)]       [ANGLEX(L)]     [ALMAX(L)]      [ALMIN(L)]      [ATMAX(L)]      [ATMIN(L)]
inp.pmaxfa = 1.;
inp.pminfa = 1.;
inp.angfac = 1.;
inp.almaxf = 1.;
inp.alminf = 1.;
inp.atmaxf = 1.;
inp.atminf = 1.;
inp.l      = (1:ne)';
inp.lreg   = zeros(ne,1)+1;


pmax_mtx_gravity_compensated_m2= zeros(size(x_ele_mtx_gravity_compensated_m)) + permeability_sand_m2 ;   %background permeability

mask_ele_mtx_silt_layer_gravity_compensated = y_ele_mtx_gravity_compensated_m  > -6  ;   % mask matrix, for element matrix 

pmax_mtx_gravity_compensated_m2 (mask_ele_mtx_silt_layer_gravity_compensated) =  permeability_silt_m2;   % silt layer permeability

pmax_mtx_m2=flip(pmax_mtx_gravity_compensated_m2);

pmax_array_m2 = pmax_mtx_m2(:);



inp.pmax   = pmax_array_m2;
inp.pmin   = pmax_array_m2;
inp.anglex = zeros(ne,1);
inp.almax  = zeros(ne,1)+0.5e-0;
inp.almin  = zeros(ne,1)+0.5e-0;
inp.atmax  = zeros(ne,1)+0.5e-0;
inp.atmin  = zeros(ne,1)+0.5e-0;

% note: for SUTRASET when node number is negative, the second input is surface area of the node
% and the third input is the thickness of the cell.

%% DATASET 17   Data for Fluid Source and Sinks

%inp.iqcp  = iqcp;
%inp.qinc  = qinc;
%inp.uinc  = uinc;


% ## DATASET 19:  Data for Specified Pressure Nodes
%###  [IPBC]                [PBC]                [UBC]

mask_nod_mtx_aquifer_boundary_gravity_compensated = and(y_nod_mtx_gravity_compensated_m<-4, x_nod_mtx_gravity_compensated_m>99.99);  % below 4 metre, greater than 200 m away from the centre
ipbc_node_idx_array                           = node_index_mtx_gravity_compensated(mask_nod_mtx_aquifer_boundary_gravity_compensated);
pbc                                           = -(y_nod_mtx_gravity_compensated_m(mask_nod_mtx_aquifer_boundary_gravity_compensated) - initial_head_aquifer_m ) *c.g * (inp.rhow0 + inp.drwdu * c_saltwater_kgPkg);



%mask_mtx_aquifer_boundary_gravity_compensated_left = and(y_nod_mtx_gravity_compensated_m<-4, x_nod_mtx_gravity_compensated_m<0.01);
%ipbc_node_idx_array_left                           = node_index_mtx_gravity_compensated(mask_mtx_aquifer_boundary_gravity_compensated_left);
%pbc_left                                      = -(y_nod_mtx_gravity_compensated_m(mask_mtx_aquifer_boundary_gravity_compensated) - initial_head_aquifer_m +0.5 ) *c.g * inp.rhow0 ;
%
%inp.ipbc = [ipbc_node_idx_array; ipbc_node_idx_array_left];
%inp.pbc  = [pbc;pbc_left];
%inp.ubc  = [zeros(size(pbc))+c_saltwater_kgPkg;zeros(size(pbc_left))+c_freshwater_kgPkg];
%inp.npbc = length(inp.pbc);



mask_nod_mtx_aquifer_boundary_gravity_compensated_top = and(y_nod_mtx_gravity_compensated_m>-0.1, x_nod_mtx_gravity_compensated_m<pond_radius_m);
ipbc_node_idx_array_top                           = node_index_mtx_gravity_compensated(mask_nod_mtx_aquifer_boundary_gravity_compensated_top);
pbc_top                                      = zeros(size(ipbc_node_idx_array_top));

%inp.ipbc = [ipbc_node_idx_array; ipbc_node_idx_array_top];
%inp.pbc  = [pbc;pbc_top];
%inp.ubc  = [zeros(size(pbc))+c_saltwater_kgPkg;zeros(size(pbc_top))+c_freshwater_kgPkg];
%inp.npbc = length(inp.pbc);



inp.ipbc = ipbc_node_idx_array;
inp.pbc  = pbc;
inp.ubc  = zeros(size(pbc))+c_saltwater_kgPkg;
inp.npbc = length(inp.pbc);


%##
%##  DATASET     22:  Ele        ment Incid      ence Data
%##    [LL]      [IIN(1)]        [IIN(2)]        [IIN(3)]        [IIN(4)]

%ne_mtx=reshape(l,ney,nex);
inp.iin1 = zeros(inp.ne,1);
inp.iin2 = zeros(inp.ne,1);
inp.iin3 = zeros(inp.ne,1);
inp.iin4 = zeros(inp.ne,1);

% DATASET 22, user need to check if the sequence of the node is in clockwise manner as suggested by the manual

%for j  = 1:nex
%    for i = 1:ney
%        inp.iin1(idx) = node_index_mtx(i,j);
%        inp.iin2(idx) = node_index_mtx(i,j+1);
%        inp.iin3(idx) = node_index_mtx(i+1,j+1);
%        inp.iin4(idx) = node_index_mtx(i+1,j);
%        idx           = idx+1;
%    end
%end


inp.iin1= iin1;
inp.iin2= iin2;
inp.iin3= iin3;
inp.iin4= iin4;

inp.export_to_file();


% setting the initial pressure as hydrostatic, in particular at the sandy aquifer, the silt layer will be overwritten
pm1_mtx_gravity_compensated_pa= - (- initial_head_aquifer_m + y_nod_mtx_gravity_compensated_m)*c.g*c.rhow_pure_water;


mask_nod_mtx_silt_layer_gravity_compensated = y_nod_mtx_gravity_compensated_m  > -6  ;   % mask matrix, for nod matrix 

%set a relatively high suction in the silt layer
pm1_mtx_gravity_compensated_pa ( mask_nod_mtx_silt_layer_gravity_compensated) = -20000;




% put the top cell as zero pressure 
pm1_mtx_gravity_compensated_pa(mask_nod_mtx_aquifer_boundary_gravity_compensated_top)=initial_pond_water_depth_m*c.rhow_pure_water*c.g;

pm1_mtx_pa=flip(pm1_mtx_gravity_compensated_pa);


% initial solute concentrtion, sandy aquifer has a saline concentration while silt has a fresh concentration
um1_mtx_gravity_compensated_kgPkg= zeros(size(pm1_mtx_gravity_compensated_pa))+ c_saltwater_kgPkg;



um1_mtx_gravity_compensated_kgPkg (mask_nod_mtx_silt_layer_gravity_compensated) = c_freshwater_kgPkg;


um1_mtx_kgPkg=flip(um1_mtx_gravity_compensated_kgPkg);





fprintf('use the original ics file')
% ics file
ics       = icsObj('PART1','read_from_file','no');
ics.tics  = 0;
%ics.cpuni = 'UNIFORM';
%ics.pm1   = -30;
ics.cpuni = 'NONUNIFORM';
ics.pm1   = pm1_mtx_pa(:);
ics.cuuni = 'NONUNIFORM';
ics.um1   = um1_mtx_kgPkg ;
ics.export_to_file();
%
%
%%pm1_mtx=reshape(b.pm1,[inp.nn1,inp.nn2]);
%%
%%b=icsObj('PART1_cs.csv','inpObj',inp);
%%pm1_mtx=reshape(b.pm1,[inp.nn1,inp.nn2]);
%%um1_mtx=reshape(b.um1,[inp.nn1,inp.nn2]);
%%inp.get_x_nod_mtx;
%%inp.get_y_nod_mtx;
%%a=figure;
%%surf(inp.x_nod_mtx,inp.y_nod_mtx,um1_mtx);
%%savefig(a,'conc.fig');
%%
%%a=figure;
%%surf(inp.x_nod_mtx,inp.y_nod_mtx,pm1_mtx);
%%savefig(a,'p.fig');
%%
%%saveas(a,'Barchart.png')
%
%
%ics.export_to_file();
