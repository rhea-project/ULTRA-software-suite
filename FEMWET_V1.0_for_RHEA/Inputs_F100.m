function [AC, X] = Inputs()

% clear all
% close all
% clc

%% load initial data from EMWET

cd Data
load AC
load Ts
cd ..

%% Wing geometry definition
% ******** FOKKER 100 Wing ***************

AC.Wing.Geom = [0         0         0         5.6000    0
                2.3400    4.6000    0.2000    3.6000    -0.98*pi/180;
                5.5000   14.0400    0.6100    1.2600    -3*pi/180];

            
AC.Wing.Airfoils.eta = [0;   4.6/14.04;    10/14.04;    1];  
     
% airfoils parallel to freestream velocity
Airfoils =  [0.1700    0.1000    0.2468    0.0982    0.1242   -0.1498   -0.1028   -0.2655   -0.1142   -0.0574
             0.1850    0.1132    0.2967    0.1278    0.2244   -0.1273   -0.0021   -0.3089   -0.0648    0.0068
             0.2041    0.1038    0.2392    0.1616    0.1941   -0.1128    0.1291   -0.3322    0.1686   -0.1491
             0.2127    0.0750    0.2441    0.1179    0.1461   -0.1304    0.0686   -0.2754    0.0533   -0.1333];


AC.Aero.Xs = 0.5;                    % sweep line
Sweep = atan((AC.Wing.Geom(end,1)+ AC.Aero.Xs*AC.Wing.Geom(end,4)-AC.Aero.Xs *AC.Wing.Geom(1,4))/AC.Wing.Geom(end,2));
AC.Wing.Sweep = Sweep;


%% Chebyshev modes for airfoil design

AC.Wing.Airfoils.nAirfoil = 8;                                          % number of sections that are defined as design variables
AC.Wing.Airfoils.yAirfoil = linspace(0,1,AC.Wing.Airfoils.nAirfoil);

AC.Wing.Airfoils.CST = interp1(AC.Wing.Airfoils.eta,Airfoils,AC.Wing.Airfoils.yAirfoil)/cos(Sweep);    % airfoils perpendicular to sweep line


% Airfoils Chebyshev perturbation
AC.Wing.Airfoils.ModeNum  = 5;  % num,ber of modes for each side
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

%% Aerodynamics

% flight conditions
H = 10675;                            % Altitude [m]
AC.Aero.Mach = 0.75;                   % Mach number

[T, a, ~, rho] = atmosisa(H);           
AC.Aero.Vinf = AC.Aero.Mach*a;        % speed
AC.Aero.rho = rho;                    % air density
AC.Aero.mu  = 1.716e-5*(273.15+110.3)/(T+110.3)*(T/273.15)^(3/2);   % Sutherland's law

%% structure

AC.St.Spars = [0         0.1000    0.6100
               0.3316    0.1600    0.6300
               1.0000    0.2200    0.5500];


AC.FEM.nNodea = AC.Q3D.Grid.DnSpan;   % node distribution due to aileron

AC.FEM.nNode  = sum(AC.FEM.nNodea);           % number of nodes
AC.FEM.nElm   = AC.FEM.nNode - 1;           % number of elements

AC.FEM.nDOF   = 6*AC.FEM.nNode;   % number of degree of freedom

AC.FEM.nSkin = 4;  % number of elements for each skin panel 
AC.FEM.nSpar = 2;  % number of elements for each spar panel

%% Weight

AC.Weight.MTOW = 43090;
AC.Weight.ZFW  = 35830;
AC.Weight.MLW  = 38780;

% AC.Weight.Wdes = 2.5*AC.Weight.MTOW*9.81;
% AC.Weight.Wdes = sqrt(AC.Weight.MTOW*AC.Weight.ZFW)*9.81;  % cruise


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


% palnform
Cr = AC.Wing.Geom(1,4);
lambda1 = AC.Wing.Geom(2,4)/AC.Wing.Geom(1,4);
lambda2 = AC.Wing.Geom(3,4)/AC.Wing.Geom(2,4);
b1 = AC.Wing.Geom(2,2);
b2 = AC.Wing.Geom(3,2)-b1;
Gamma1 = atan(AC.Wing.Geom(2,1)/b1);
Gamma2 = atan((AC.Wing.Geom(3,1)-AC.Wing.Geom(2,1))/b2);
epsilon1 = AC.Wing.Geom(2,end)*pi/180;
epsilon2 = AC.Wing.Geom(3,end)*pi/180;

P = [Cr; lambda1; lambda2; b1; b2; Gamma1; Gamma2; epsilon1; epsilon2];

% Design Vector
X = [T; reshape(AC.Wing.Airfoils.Chebyshev,size(AC.Wing.Airfoils.Chebyshev,1)*size(AC.Wing.Airfoils.Chebyshev,2),1);P];

% **************** TO Be Implemented  ****************
% 1- effect of twist on structure
% 2- effect of dihedral
% 3- fuel loads
% 4- Left aileron modeling
% ****************************************************
