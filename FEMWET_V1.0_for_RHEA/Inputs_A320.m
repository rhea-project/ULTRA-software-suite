function [AC, T, G, P] = Inputs_A320()

% clear all
% close all
% clc

%% load initial data from EMWET

cd Data
load T_A320
cd ..

%% Wing geometry definition
% ******** FOKKER 100 Wing ***************

AC.Wing.Geom = [0        0         0        7.0518   0
                3.3006   6.3403    0.5547   3.7584   -2.5
                8.8306   16.9635   1.4841   1.4958   -2.5];

AC.Wing.kink = 6.3403/16.9635;
            
AC.Wing.Airfoils.eta = [0;   0.33;    0.66;    1];  
     
% airfoils parallel to freestream velocity
Airfoils =  [0.2021    0.1175    0.2342    0.2067    0.2841   -0.1265   -0.1101   -0.3012   -0.1164    0.1813
             0.1849    0.1072    0.2216    0.1780    0.2817   -0.1529    0.0026   -0.4241    0.1435   -0.0779
             0.1635    0.1116    0.1996    0.1664    0.2751   -0.1617    0.0170   -0.4009    0.1635   -0.0979
             0.1542    0.1051    0.1955    0.1603    0.2673   -0.1410   -0.0371   -0.2335   -0.0791    0.1702];


AC.Aero.Xs = 0.5;                    % sweep line
Sweep = atan((AC.Wing.Geom(end,1)+ AC.Aero.Xs*AC.Wing.Geom(end,4)-AC.Aero.Xs *AC.Wing.Geom(1,4))/AC.Wing.Geom(end,2));
AC.Wing.Sweep = Sweep;


%% Chebyshev modes for airfoil design

AC.Wing.Airfoils.nAirfoil = 8;                                          % number of sections that are defined as design variables
AC.Wing.Airfoils.yAirfoil = linspace(0,1,AC.Wing.Airfoils.nAirfoil);

AC.Wing.Airfoils.CST = interp1(AC.Wing.Airfoils.eta,Airfoils,AC.Wing.Airfoils.yAirfoil)/cos(Sweep);    % airfoils perpendicular to sweep line


% Airfoils Chebyshev perturbation
AC.Wing.Airfoils.ModeNum  = 10;  % number of modes for each side
AC.Wing.Airfoils.Chebyshev = zeros(AC.Wing.Airfoils.nAirfoil,2*AC.Wing.Airfoils.ModeNum);


%% Aileron geometry

%                  eta_i   eta_f    cf/c   delta_f (downward) 
AC.Aileron.Geom = [ 0.75    0.95    0.25    1e-6 ];     % each row for one flap



%% Q3D setting

% grids
AC.Q3D.Grid.DnSpan = [7 3 1];          % distribution of spanwise grids
AC.Q3D.Grid.nSpan = sum(AC.Q3D.Grid.DnSpan) - 1;    % number of spanwise grids
% AC.Q3D.Grid.nPanel = AC.Q3D.Grid.nSpan;

AC.Q3D.Grid.DnChord = [4 2];            % Distribution of chordwise grids
AC.Q3D.Grid.nChord = sum(AC.Q3D.Grid.DnChord) - 1;  % number of chordwise grids

AC.Q3D.N2D = 8;                    % number of 2D sections to be analysed

%% Calibration factors

Kcdi = 1;
Kcdp = 1;
Kcdf = 1;

Kcd = [Kcdi Kcdp Kcdf];

cd FEMWET
save('Kcd.mat','Kcd')
cd ..


%% structure and Materials

AC.St.Spars = [0         0.1000    0.6100
               0.3738    0.1600    0.6300
               1.0000    0.2200    0.5500];

AC.St.FuelTank = [0.05 0.85];   % start and end of fuel tank           

AC.FEM.nNodea = AC.Q3D.Grid.DnSpan;   % node distribution due to aileron

AC.FEM.nNode  = sum(AC.FEM.nNodea);           % number of nodes
AC.FEM.nElm   = AC.FEM.nNode - 1;           % number of elements

AC.FEM.nDOF   = 6*AC.FEM.nNode;   % number of degree of freedom

AC.FEM.nSkin = 4;  % number of elements for each skin panel 
AC.FEM.nSpar = 2;  % number of elements for each spar panel

cd Data
AC.St.Mat_u  = importdata('Al_7075_T6.mat');
AC.St.Mat_l  = importdata('Al_2024_T3.mat');
AC.St.Mat_fs = importdata('Al_7075_T6.mat');
AC.St.Mat_rs = importdata('Al_7075_T6.mat');
cd ..

AC.St.F =  0.9600;  % stiffened panel efficiency
AC.St.RP = 0.5867;  % rib pitch

%% Weight

AC.Weight.MTOW = 73500;
AC.Weight.FW   = 17940;
AC.Weight.ZFW  = 73500-17940;
AC.Weight.MLW  = 0.8*73500;

% AC.Weight.Wdes = 2.5*AC.Weight.MTOW*9.81;
% AC.Weight.Wdes = sqrt(AC.Weight.MTOW*AC.Weight.ZFW)*9.81;  % cruise

%% Load cases

% Case1 2.5g pull up
AC.Load{1}.H = 7500; 
AC.Load{1}.Mach = 0.89;
AC.Load{1}.n = 2.5;
    
% Case 2 2.5g pull up
AC.Load{2}.H = 0; 
AC.Load{2}. Mach = 0.58;
AC.Load{2}.n = 2.5;
    
% Case 3 -1g push down
AC.Load{3}.H = 7500; 
AC.Load{3}.Mach = 0.89;
AC.Load{3}.n = -1;

% Case 4 Fatigue (1.3g gust)
AC.Load{4}.H = 7500; 
AC.Load{4}.Mach = 0.89;
AC.Load{4}.n = 1.3;
    
% Case 5 Roll
AC.Load{5}.H = 4000; 
AC.Load{5}.Mach = 0.83;
AC.Load{5}.n = 1;
    
% Case 6 cruise
AC.Load{6}.H = 11000; 
AC.Load{6}.Mach = 0.78;  
AC.Load{6}.n = 1;


%% Design variables

% structure
nT = 10;
yT = linspace(0,1,nT);
T  = [interp1(Ts.yT,Ts.T(1:10),yT) interp1(Ts.yT,Ts.T(11:20),yT) interp1(Ts.yT,Ts.T(21:30),yT)  interp1(Ts.yT,Ts.T(31:40),yT)]';

for i=1:length(T)
    if T(i) <0.002
        T(i) = 0.002;
    end
end

AC.Structure.nT = nT;
AC.Structure.yT = yT;

% airfoils

G = reshape(AC.Wing.Airfoils.Chebyshev,size(AC.Wing.Airfoils.Chebyshev,1)*size(AC.Wing.Airfoils.Chebyshev,2),1);

% palnform
Cr = AC.Wing.Geom(1,4);
lambda = AC.Wing.Geom(3,4)/AC.Wing.Geom(1,4);
b = AC.Wing.Geom(3,2);
Gamma = atan(AC.Wing.Geom(3,1)/b);
epsilon1 = AC.Wing.Geom(2,5)*pi/180;
epsilon2 = AC.Wing.Geom(3,5)*pi/180;

 P = [Cr; lambda; b; Gamma; epsilon1; epsilon2];
% P = [Cr; lambda; b; Gamma];

% % jig shape
% 
% nJ = 10;
% yJ = linspace(0,1,nJ);
% Eps = interp1(AC.Wing.Geom(:,2),AC.Wing.Geom(:,end),yJ*b)*pi/180;
% Z = interp1(AC.Wing.Geom(:,2),AC.Wing.Geom(:,3),yJ*b);
% 
% P = [P; Z'; Eps'];
% 
% AC.Jig.nJ = nJ;
% AC.Jig.yJ = yJ;


% **************** TO Be Implemented  ****************
% 1- effect of twist on structure
% 2- effect of dihedral
% 3- fuel loads
% 4- Left aileron modeling
% ****************************************************
