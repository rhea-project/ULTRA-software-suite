clear all
close all
clc

warning off

%%

if exist('gradientinit')==0
    cd AD10
    Start;
    cd ..
end


% cd AD10
% startup;
% cd ..


%% Inputs

[AC, T, G, P] = Inputs_A320;


%% defining the design variables


DV = 7;    % 1  only structure variables                                    (aeroelastic optimization)
           % 2  only airfoil shape variables                                (pure aerodynamic shape optimization)
           % 3  only planform variables                                     (pure aerodynamic shape optimization)
           % 4  airfoil and planform variables                              (pure aerodynamic shape optimization)
           % 5  structure and airfoil variables                             (Aerostructural optimization)
           % 6  structure and planform variables                            (Aerostructural optimization)
           % 7  structure, airfoil and planform variables                   (Aerostructural optimization)
           

%%

cd FEMWET
tic
[Wwing, Failure, U, eta_a, CD, WL, Vf, dWw_dX, Dg_Dx,dU_dX, dMa_dX, dCD_dX, dWL_dX, dVf_dX, CDi, CDp, CDf, dM_da,dM_da0]=FEMWET(AC,T,G,P,6,1,1,0,DV,AC.Weight.FW,AC.Weight.MTOW);
t = toc
cd ..
