%clear all
%close all
% project name
name='Column';

% read input files
%inp=i0('desaline');
fil =  readFIL;
% inp  = inpObj(fil.basename,'block_reading','yes');

%%need to replace REK\n     Y velocity     Z velocity in ele to Y velocity     Z velocity     REK
nod  = readNOD(fil.basename);
ele  = readELE(fil.basename);
bcop = readBCOP(fil.basename);
bcof = readBCOF(fil.basename);

% index for p c and s
x_idx  = strcmp(nod(1).label,'X');
y_idx  = strcmp(nod(1).label,'Y');
z_idx  = strcmp(nod(1).label,'Z');

p_idx  = strcmp(nod(1).label,'Pressure');
c_idx  = strcmp(nod(1).label,'Concentration');
s_idx  = strcmp(nod(1).label,'Saturation');
sm_idx = strcmp(nod(1).label,'SM');
vz_idx = strcmp(ele(1).label,'Z velocity');
vy_idx = strcmp(ele(1).label,'Y velocity');
vx_idx = strcmp(ele(1).label,'X velocity');
rek_idx = strcmp(ele(1).label,'REK');
zele_idx = strcmp(ele(1).label,'Z origin');

%read swcc data;
swcc=swccObj('ewatering_swcc_parameters.dat');
swcc.calc_swcc('type','fayer1995');
swcc.psim=-[0:0.01:5,5:0.1:10,11:1:50,51:10:100,200:10:1000,1000:100:10000,11000:1000:50000]';
swcc.calc_swcc('type','fayer1995');
swcc.calc_relativek('type','mualem1976','tortuosity',0.5);