function [K2,F2,M2] = Element_int(Gvar,AC,Load,nNodea,Conc,ele,nT,yT,nAirfoil,yAirfoil,nC,n,nTs,nCn,DV,DP)

X = Gvar(1:nTs);
Xea = Gvar(nTs+1);
Yea = Gvar(nTs+2);
cNode = Gvar(nTs+3:nTs+2+nCn);
Vf = Gvar(nTs+3+nCn);
Wf = Gvar(nTs+4+nCn);

ModeNum  = AC.Wing.Airfoils.ModeNum;
nNode = sum(nNodea);

%%


if DV == 1
    T  = [X(1:end-2); DP; X(end-1:end)];
elseif DV == 2
    T  = [DP(1:4*nT); X(1:end-2); DP(4*nT+1:end); X(end-1:end)];
elseif DV == 3 
    T  = [DP; X];
elseif DV ==4
    T  = [DP; X];
elseif DV ==5
    T  = [X(1:end-2); DP; X(end-1:end)];
elseif DV ==6    
    T  = [X(1:4*nT); DP; X(4*nT+1:end)];
elseif DV ==7
    T  = X;
end

%%

CST      = AC.Wing.Airfoils.CST;

Cheby   = T(4*nT+1:4*nT+2*ModeNum*nAirfoil);
Cheby   = reshape(Cheby,nAirfoil,2*ModeNum);


Tu = T(1:nT); Tl = T(nT+1:2*nT); Tfs = T(2*nT+1:3*nT); Trs = T(3*nT+1:4*nT);

% CSTu = T(4*nT+1:4*nT+5*nAirfoil);
% CSTl = T(4*nT+5*nAirfoil+1:4*nT+2*5*nAirfoil);

    
Cr     = T(4*nT+2*ModeNum*nAirfoil+1);
lambda = T(4*nT+2*ModeNum*nAirfoil+2);
b      = T(4*nT+2*ModeNum*nAirfoil+3);
Gamma  = T(4*nT+2*ModeNum*nAirfoil+4);

b1 = AC.Wing.kink*b;
b2 = b - b1;

%       x                      y          z                   chord
Geom = [0                      0          AC.Wing.Geom(1,3)   Cr; 
        b1*tan(Gamma)          b1         AC.Wing.Geom(2,3)   Cr-AC.Wing.kink*b*tan(Gamma) ;
        b*tan(Gamma)           b          AC.Wing.Geom(3,3)   lambda*Cr];
  
    
    
cNode = [cNode(1:nNode) cNode(nNode+1:2*nNode) cNode(2*nNode+1:3*nNode)];

sweepEA = atan((cNode(end,1)-cNode(1,1))/(cNode(end,2)-cNode(1,2)));


%% element geometry

C1(1) = cNode(Conc(ele,1),2);   % coordinate of element start
C1(2) = cNode(Conc(ele,1),3);
C1(3) = cNode(Conc(ele,1),1);

C2(1) = cNode(Conc(ele,2),2);   % coordinate of element end
C2(2) = cNode(Conc(ele,2),3);
C2(3) = cNode(Conc(ele,2),1);

% L     = norm(C2-C1);          % Element length 
L = sqrt(sum((C2-C1).^2));

Coy   = [0 1 0];
Coz   = [0 0 1];


lox = (C2(1)-C1(1))/L;         % Find direction cosines and lambda matrix.
mox = (C2(2)-C1(2))/L;
nox = (C2(3)-C1(3))/L;

lam = [lox mox nox; Coy; Coz];
zerom = zeros(3,3);
Tr = [ lam zerom zerom zerom          % Transformation matrix
  zerom lam zerom zerom
  zerom zerom lam zerom
  zerom zerom zerom lam];

%% Defining variables and constants 

k = 5/6;

eta_ai = AC.Aileron.Geom(1); 
eta_af = AC.Aileron.Geom(2);
yN1 = linspace(0,eta_ai,nNodea(1));
yN2 = linspace(eta_ai,eta_af,nNodea(2)+1);
yN3 = linspace(eta_af,1,nNodea(3)+1);
yNode = [yN1 yN2(2:end) yN3(2:end)];

% yNode = linspace(0,1,nNode);
eta_i = yNode(Conc(ele,1));  % start of the element
eta_f = yNode(Conc(ele,2));  % end of the element

%% CST as design variables
CST_i = zeros(1,size(CST,2));
for i=1:size(CST,1)-1
    CST_i = CST_i + (meshgrid(CST(i,:),1)+ (meshgrid(CST(i+1,:),1)-meshgrid(CST(i,:),1)).*...
        meshgrid((eta_i-yAirfoil(i))/((yAirfoil(i+1)-yAirfoil(i))),1:size(CST,2))').*meshgrid(heaviside(heaviside(eta_i-(yAirfoil(i)))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1))-eta_i)-0.6),1:size(CST,2))';
end
if eta_i == yAirfoil(end)
    CST_i(end,:) = CST(end,:);
end

CST_f = zeros(1,size(CST,2));
for i=1:size(CST,1)-1
    CST_f = CST_f + (meshgrid(CST(i,:),1)+ (meshgrid(CST(i+1,:),1)-meshgrid(CST(i,:),1)).*...
        meshgrid((eta_f-yAirfoil(i))/((yAirfoil(i+1)-yAirfoil(i))),1:size(CST,2))').*meshgrid(heaviside(heaviside(eta_f-(yAirfoil(i)))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1))-eta_f)-0.6),1:size(CST,2))';
end
if eta_f == yAirfoil(end)
    CST_f(end,:) = CST(end,:);
end

Cheby_i = zeros(1,size(Cheby,2));
for i=1:size(CST,1)-1
    Cheby_i = Cheby_i + (meshgrid(Cheby(i,:),1)+ (meshgrid(Cheby(i+1,:),1)-meshgrid(Cheby(i,:),1)).*...
        meshgrid((eta_i-yAirfoil(i))/((yAirfoil(i+1)-yAirfoil(i))),1:size(Cheby,2))').*meshgrid(heaviside(heaviside(eta_i-(yAirfoil(i)))-0.1),1:size(Cheby,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1))-eta_i)-0.6),1:size(Cheby,2))';
end
if eta_i == yAirfoil(end)
    Cheby_i(end,:) = Cheby(end,:);
end

Cheby_f = zeros(1,size(Cheby,2));
for i=1:size(CST,1)-1
    Cheby_f = Cheby_f + (meshgrid(Cheby(i,:),1)+ (meshgrid(Cheby(i+1,:),1)-meshgrid(Cheby(i,:),1)).*...
        meshgrid((eta_f-yAirfoil(i))/((yAirfoil(i+1)-yAirfoil(i))),1:size(Cheby,2))').*meshgrid(heaviside(heaviside(eta_f-(yAirfoil(i)))-0.1),1:size(Cheby,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1))-eta_f)-0.6),1:size(Cheby,2))';
end
if eta_f == yAirfoil(end)
    Cheby_f(end,:) = Cheby(end,:);
end

%% Wing geometry at the start and end of the element

Chord_i = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*cNode(Conc(ele,1),2)/b1).*heaviside(heaviside(b1-cNode(Conc(ele,1),2))-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(cNode(Conc(ele,1),2)-b1)/b2).*heaviside(heaviside(cNode(Conc(ele,1),2)-b1)-0.6);

Chord_i = Chord_i*cos(sweepEA);  % section perpendicular to the sweep line

Chord_f = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*cNode(Conc(ele,2),2)/b1).*heaviside(heaviside(b1-cNode(Conc(ele,2),2))-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(cNode(Conc(ele,2),2)-b1)/b2).*heaviside(heaviside(cNode(Conc(ele,2),2)-b1)-0.6);

Chord_f = Chord_f*cos(sweepEA);  % section perpendicular to the sweep line

for i=1:length(AC.St.Spars(:,1));
    if eta_i == AC.St.Spars(i,1)
        fs_i = AC.St.Spars(i,2);
        rs_i = AC.St.Spars(i,3);
                break;
    elseif eta_i > AC.St.Spars(i,1) && eta_i < AC.St.Spars(i+1,1)
        fs_i = AC.St.Spars(i,2) + (AC.St.Spars(i+1,2)-AC.St.Spars(i,2))*(eta_i-AC.St.Spars(i,1))/(AC.St.Spars(i+1,1)-AC.St.Spars(i,1));
        rs_i = AC.St.Spars(i,3) + (AC.St.Spars(i+1,3)-AC.St.Spars(i,3))*(eta_i-AC.St.Spars(i,1))/(AC.St.Spars(i+1,1)-AC.St.Spars(i,1));
                break;
    end
end

for i=1:length(AC.St.Spars(:,1));
    if eta_f == AC.St.Spars(i,1)
        fs_f = AC.St.Spars(i,2);
        rs_f = AC.St.Spars(i,3);
                break;
    elseif eta_f > AC.St.Spars(i,1) && eta_f < AC.St.Spars(i+1,1)
        fs_f = AC.St.Spars(i,2) + (AC.St.Spars(i+1,2)-AC.St.Spars(i,2))*(eta_f-AC.St.Spars(i,1))/(AC.St.Spars(i+1,1)-AC.St.Spars(i,1));
        rs_f = AC.St.Spars(i,3) + (AC.St.Spars(i+1,3)-AC.St.Spars(i,3))*(eta_f-AC.St.Spars(i,1))/(AC.St.Spars(i+1,1)-AC.St.Spars(i,1));
                break;
    end
end


for i=1:nT
    if eta_i == yT(i)
        tu_i = Tu(i);
        tl_i = Tl(i);
        tfs_i = Tfs(i);
        trs_i = Trs(i);
                break;
    elseif eta_i > yT(i) && eta_i < yT(i+1)
        tu_i =  Tu(i) + (Tu(i+1)-Tu(i))*(eta_i-yT(i))/(yT(i+1)-yT(i));
        tl_i =  Tl(i) + (Tl(i+1)-Tl(i))*(eta_i-yT(i))/(yT(i+1)-yT(i));
        tfs_i =  Tfs(i) + (Tfs(i+1)-Tfs(i))*(eta_i-yT(i))/(yT(i+1)-yT(i));
        trs_i =  Trs(i) + (Trs(i+1)-Trs(i))*(eta_i-yT(i))/(yT(i+1)-yT(i));
                break;
    end
end

for i=1:nT
    if eta_f == yT(i)
        tu_f = Tu(i);
        tl_f = Tl(i);
        tfs_f = Tfs(i);
        trs_f = Trs(i);
                break;
    elseif eta_f > yT(i) && eta_f < yT(i+1)
        tu_f =  Tu(i) + (Tu(i+1)-Tu(i))*(eta_f-yT(i))/(yT(i+1)-yT(i));
        tl_f =  Tl(i) + (Tl(i+1)-Tl(i))*(eta_f-yT(i))/(yT(i+1)-yT(i));
        tfs_f =  Tfs(i) + (Tfs(i+1)-Tfs(i))*(eta_f-yT(i))/(yT(i+1)-yT(i));
        trs_f =  Trs(i) + (Trs(i+1)-Trs(i))*(eta_f-yT(i))/(yT(i+1)-yT(i));
                break;
    end
end

%% Structural properties

% Modified for new coordinate system (x,y,z) <---- (y, z, x)   
% Iy = Iz (Iy in EMWET definition) and  Iz = Ix;


[~,~,~,~, Aa, Aza, Aya, Izza, Iyya, Iyza, Ja, EAa, EAza, EAya,  EIzza, EIyya, EIyza, GAa, GAza, GAya, GJa, Va] ... 
    = Section_Properties(CST_i,Cheby_i,fs_i,rs_i,Chord_i,AC.St.Mat_u,AC.St.Mat_l,AC.St.Mat_fs,AC.St.Mat_rs,tu_i,tl_i,tfs_i,trs_i,Xea,Yea);

[~,~,~,~, Ab, Azb, Ayb, Izzb, Iyyb, Iyzb, Jb, EAb, EAzb, EAyb,  EIzzb, EIyyb, EIyzb, GAb, GAzb, GAyb, GJb, Vb] ... 
    = Section_Properties(CST_f,Cheby_f,fs_f,rs_f,Chord_f,AC.St.Mat_u,AC.St.Mat_l,AC.St.Mat_fs,AC.St.Mat_rs,tu_f,tl_f,tfs_f,trs_f,Xea,Yea);


 %% wing and fuel weight relief
 
r = 2800*2;  % material density (multiplied by 2 for secondary weights)
Paw = -n*Aa*r*9.81;
Pbw = -n*Ab*r*9.81;

Vfe = 0.93*(Va+Vb)/2*L;

if eta_i >=AC.St.FuelTank(1) && eta_f <=AC.St.FuelTank(2)
    Wfe = Wf*Vfe/Vf;
else
    Wfe = 0;
end

Paf = -n*2*Wfe/L/(1+Vb/Va);
Pbf = -n*2*Wfe/L/(1+Va/Vb);


%% External loads

Gamma = Load(ele:nNode-1:end);

Load2 = [Gamma(1); Gamma(2:end)-Gamma(1:end-1)];
Lift = sum(Load2);

etaC = linspace(0,1,nC+1);
etaCc = etaC(1:end-1) + 0.75*(etaC(2:end)-etaC(1:end-1));

Pa = Lift + Paw + Paf;  % assuming constant lift
Pb = Lift + Pbw + Pbf;

Ta = sum(Load2.*(Xea-etaCc)'*Chord_i)/cos(sweepEA);
Tb = sum(Load2.*(Xea-etaCc)'*Chord_f)/cos(sweepEA);

%% Load interpolation

% Gamma = Load(nNode:end)-Load(1:end-nNode+1);
% Gamma = [Load(1:nS-1); Gamma];
% 
% Lift = Gamma(end-nNode+2:end);

%% Stiffness matrix
% instotrpic



K1 = [    (EAa + EAb)/(2*L)
             -(EAza - EAzb)/L^2
              (EAya - EAyb)/L^2
                              0
                        -EAya/L
                        -EAza/L
             -(EAa + EAb)/(2*L)
              (EAza - EAzb)/L^2
             -(EAya - EAyb)/L^2
                              0
                         EAyb/L
                         EAzb/L];
         
K2 = [      -(EAza - EAzb)/L^2
           (6*(EIzza + EIzzb))/L^3
          -(6*(EIyza + EIyzb))/L^3
                                 0
         (2*(2*EIyza + EIyzb))/L^2
         (2*(2*EIzza + EIzzb))/L^2
                 (EAza - EAzb)/L^2
          -(6*(EIzza + EIzzb))/L^3
           (6*(EIyza + EIyzb))/L^3
                                 0
         (2*(EIyza + 2*EIyzb))/L^2
         (2*(EIzza + 2*EIzzb))/L^2];

K3 = [                                  (EAya - EAyb)/L^2
                                     -(6*(EIyza + EIyzb))/L^3
  (((12*GAa*k)/5 + (12*GAb*k)/5)*L^2 + 6*EIyya + 6*EIyyb)/L^3
                                              (GAza + GAzb)/L
                                   -(2*(2*EIyya + EIyyb))/L^2
                                   -(2*(2*EIyza + EIyzb))/L^2
                                           -(EAya - EAyb)/L^2
                                      (6*(EIyza + EIyzb))/L^3
 -(((12*GAa*k)/5 + (12*GAb*k)/5)*L^2 + 6*EIyya + 6*EIyyb)/L^3
                                             -(GAza + GAzb)/L
                      - (2*GAa*k)/5 - (2*EIyya + 4*EIyyb)/L^2
                                   -(2*(EIyza + 2*EIyzb))/L^2];
                               
K4 = [        0
                  0
    (GAza + GAzb)/L
  (GJa + GJb)/(2*L)
                  0
                  0
                  0
                  0
   -(GAza + GAzb)/L
 -(GJa + GJb)/(2*L)
    GAzb/6 - GAza/6
                  0];
              
K5 = [                               -EAya/L
                       (2*(2*EIyza + EIyzb))/L^2
         - (2*GAb*k)/5 - (4*EIyya + 2*EIyyb)/L^2
                                 GAza/6 - GAzb/6
                             (3*EIyya + EIyyb)/L
                             (3*EIyza + EIyzb)/L
                                          EAya/L
                      -(2*(2*EIyza + EIyzb))/L^2
           (2*GAb*k)/5 + (4*EIyya + 2*EIyyb)/L^2
                                 GAzb/6 - GAza/6
 (EIyya + EIyyb)/L - L*((GAa*k)/15 + (GAb*k)/15)
                               (EIyza + EIyzb)/L];
                           

K6 = [                  -EAza/L
          (2*(2*EIzza + EIzzb))/L^2
         -(2*(2*EIyza + EIyzb))/L^2
                                  0
                (3*EIyza + EIyzb)/L
                (3*EIzza + EIzzb)/L
                             EAza/L
         -(2*(2*EIzza + EIzzb))/L^2
          (2*(2*EIyza + EIyzb))/L^2
                                  0
                  (EIyza + EIyzb)/L
                  (EIzza + EIzzb)/L];
              
          
K7 = [   -(EAa + EAb)/(2*L)
              (EAza - EAzb)/L^2
             -(EAya - EAyb)/L^2
                              0
                         EAya/L
                         EAza/L
              (EAa + EAb)/(2*L)
             -(EAza - EAzb)/L^2
              (EAya - EAyb)/L^2
                              0
                        -EAyb/L
                        -EAzb/L];
                    
        
K8 = [            (EAza - EAzb)/L^2
               -(6*(EIzza + EIzzb))/L^3
                (6*(EIyza + EIyzb))/L^3
                                      0
             -(2*(2*EIyza + EIyzb))/L^2
             -(2*(2*EIzza + EIzzb))/L^2
                     -(EAza - EAzb)/L^2
                (6*(EIzza + EIzzb))/L^3
               -(6*(EIyza + EIyzb))/L^3
                                      0
             -(2*(EIyza + 2*EIyzb))/L^2
             -(2*(EIzza + 2*EIzzb))/L^2];


K9 = [                                 -(EAya - EAyb)/L^2
                                      (6*(EIyza + EIyzb))/L^3
 -(((12*GAa*k)/5 + (12*GAb*k)/5)*L^2 + 6*EIyya + 6*EIyyb)/L^3
                                             -(GAza + GAzb)/L
                                    (2*(2*EIyya + EIyyb))/L^2
                                    (2*(2*EIyza + EIyzb))/L^2
                                            (EAya - EAyb)/L^2
                                     -(6*(EIyza + EIyzb))/L^3
  (((12*GAa*k)/5 + (12*GAb*k)/5)*L^2 + 6*EIyya + 6*EIyyb)/L^3
                                              (GAza + GAzb)/L
                        (2*GAa*k)/5 + (2*EIyya + 4*EIyyb)/L^2
                                    (2*(EIyza + 2*EIyzb))/L^2];
  
                                
K10 = [               0
                          0
           -(GAza + GAzb)/L
         -(GJa + GJb)/(2*L)
                          0
                          0
                          0
                          0
            (GAza + GAzb)/L
          (GJa + GJb)/(2*L)
            GAza/6 - GAzb/6
                          0];
                      
                      
K11  = [                                       EAyb/L
                                (2*(EIyza + 2*EIyzb))/L^2
                  - (2*GAa*k)/5 - (2*EIyya + 4*EIyyb)/L^2
                                          GAzb/6 - GAza/6
                                        (EIyya + EIyyb)/L
                                        (EIyza + EIyzb)/L
                                                  -EAyb/L
                               -(2*(EIyza + 2*EIyzb))/L^2
                    (2*GAa*k)/5 + (2*EIyya + 4*EIyyb)/L^2
                                          GAza/6 - GAzb/6
     L*((2*GAa*k)/15 + (2*GAb*k)/5) + (EIyya + 3*EIyyb)/L
                                      (EIyza + 3*EIyzb)/L];
                              
K12 = [                  EAzb/L
          (2*(EIzza + 2*EIzzb))/L^2
         -(2*(EIyza + 2*EIyzb))/L^2
                                  0
                  (EIyza + EIyzb)/L
                  (EIzza + EIzzb)/L
                            -EAzb/L
         -(2*(EIzza + 2*EIzzb))/L^2
          (2*(EIyza + 2*EIyzb))/L^2
                                  0
                (EIyza + 3*EIyzb)/L
                (EIzza + 3*EIzzb)/L];
    

K = [K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12];

%% Load matrix 

F=[                    0
    (L*(7*Pa + 3*Pb))/20
                       0
       (L*(2*Ta + Tb))/6
                       0
  (L^2*(3*Pa + 2*Pb))/60
                       0
    (L*(3*Pa + 7*Pb))/20
                       0
       (L*(Ta + 2*Tb))/6
                       0
 -(L^2*(2*Pa + 3*Pb))/60];


%% Mass Matrix

% M1 = [    (L*r*(3*Aa + Ab))/12
%               (r*(3*Aza + 2*Azb))/10
%              -(r*(3*Aya + 2*Ayb))/10
%                                    0
%              -(L*r*(6*Aya - Ayb))/60
%              -(L*r*(6*Aza - Azb))/60
%                   (L*r*(Aa + Ab))/12
%              -(r*(3*Aza + 2*Azb))/10
%               (r*(3*Aya + 2*Ayb))/10
%                                    0
%               (L*r*(4*Aya + Ayb))/60
%               (L*r*(4*Aza + Azb))/60];
% 
% M2 = [                             (r*(3*Aza + 2*Azb))/10
%       (r*(21*Izza + 21*Izzb))/(35*L) + (L*r*(10*Aa + 3*Ab))/35
%                                     -(3*r*(Iyza + Iyzb))/(5*L)
%                                     -(L*r*(16*Aya + 5*Ayb))/60
%                                                    (Iyzb*r)/10
%                       (r*(42*Izzb + 15*Aa*L^2 + 7*Ab*L^2))/420
%                                         (r*(2*Aza + 3*Azb))/10
%  (3*L*r*(3*Aa + 3*Ab))/140 - (3*r*(28*Izza + 28*Izzb))/(140*L)
%                                      (3*r*(Iyza + Iyzb))/(5*L)
%                                      -(L*r*(5*Aya + 4*Ayb))/60
%                                                    (Iyza*r)/10
%                       -(r*(7*Aa*L^2 - 42*Izza + 6*Ab*L^2))/420];
%                   
%                   
% M3 = [                            -(r*(3*Aya + 2*Ayb))/10
%                                     -(3*r*(Iyza + Iyzb))/(5*L)
%       (r*(21*Iyya + 21*Iyyb))/(35*L) + (L*r*(10*Aa + 3*Ab))/35
%                                      (L*r*(16*Aza + 5*Azb))/60
%                      -(r*(42*Iyyb + 15*Aa*L^2 + 7*Ab*L^2))/420
%                                                   -(Iyzb*r)/10
%                                        -(r*(2*Aya + 3*Ayb))/10
%                                      (3*r*(Iyza + Iyzb))/(5*L)
%  (3*L*r*(3*Aa + 3*Ab))/140 - (3*r*(28*Iyya + 28*Iyyb))/(140*L)
%                                       (L*r*(5*Aza + 4*Azb))/60
%                        (r*(7*Aa*L^2 - 42*Iyya + 6*Ab*L^2))/420
%                                                   -(Iyza*r)/10];
%                                               
% M4= [                         0
%          -(L*r*(16*Aya + 5*Ayb))/60
%           (L*r*(16*Aza + 5*Azb))/60
%                (L*r*(3*Ja + Jb))/12
%           -(L^2*r*(2*Aza + Azb))/60
%           -(L^2*r*(2*Aya + Ayb))/60
%                                   0
%           -(L*r*(4*Aya + 5*Ayb))/60
%            (L*r*(4*Aza + 5*Azb))/60
%                  (L*r*(Ja + Jb))/12
%              (L^2*r*(Aza + Azb))/60
%              (L^2*r*(Aya + Ayb))/60];
%          
% M5 = [                    -(L*r*(6*Aya - Ayb))/60
%                                           (Iyzb*r)/10
%             -(r*(42*Iyyb + 15*Aa*L^2 + 7*Ab*L^2))/420
%                             -(L^2*r*(2*Aza + Azb))/60
%   (L*r*(84*Iyya + 28*Iyyb + 5*Aa*L^2 + 3*Ab*L^2))/840
%                              (L*r*(3*Iyza + Iyzb))/30
%                                (L*r*(Aya + 4*Ayb))/60
%                                          -(Iyzb*r)/10
%              -(r*(6*Aa*L^2 - 42*Iyyb + 7*Ab*L^2))/420
%                               -(L^2*r*(Aza + Azb))/60
%  -(L*r*(14*Iyya + 14*Iyyb + 3*Aa*L^2 + 3*Ab*L^2))/840
%                               -(L*r*(Iyza + Iyzb))/60];
%                           
%                           
% M6 = [                    -(L*r*(6*Aza - Azb))/60
%              (r*(42*Izzb + 15*Aa*L^2 + 7*Ab*L^2))/420
%                                          -(Iyzb*r)/10
%                             -(L^2*r*(2*Aya + Ayb))/60
%                              (L*r*(3*Iyza + Iyzb))/30
%   (L*r*(84*Izza + 28*Izzb + 5*Aa*L^2 + 3*Ab*L^2))/840
%                                (L*r*(Aza + 4*Azb))/60
%               (r*(6*Aa*L^2 - 42*Izzb + 7*Ab*L^2))/420
%                                           (Iyzb*r)/10
%                               -(L^2*r*(Aya + Ayb))/60
%                               -(L*r*(Iyza + Iyzb))/60
%  -(L*r*(14*Izza + 14*Izzb + 3*Aa*L^2 + 3*Ab*L^2))/840];
% 
% 
% M7 = [        (L*r*(Aa + Ab))/12
%               (r*(2*Aza + 3*Azb))/10
%              -(r*(2*Aya + 3*Ayb))/10
%                                    0
%               (L*r*(Aya + 4*Ayb))/60
%               (L*r*(Aza + 4*Azb))/60
%                 (L*r*(Aa + 3*Ab))/12
%              -(r*(2*Aza + 3*Azb))/10
%               (r*(2*Aya + 3*Ayb))/10
%                                    0
%               (L*r*(Aya - 6*Ayb))/60
%               (L*r*(Aza - 6*Azb))/60];
% 
% 
% M8 = [                                       -(r*(3*Aza + 2*Azb))/10
%              (3*L*r*(3*Aa + 3*Ab))/140 - (3*r*(28*Izza + 28*Izzb))/(140*L)
%                                                  (3*r*(Iyza + Iyzb))/(5*L)
%                                                  -(L*r*(4*Aya + 5*Ayb))/60
%                                                               -(Iyzb*r)/10
%                                    (r*(6*Aa*L^2 - 42*Izzb + 7*Ab*L^2))/420
%                                                    -(r*(2*Aza + 3*Azb))/10
%                   (r*(21*Izza + 21*Izzb))/(35*L) + (L*r*(3*Aa + 10*Ab))/35
%                                                 -(3*r*(Iyza + Iyzb))/(5*L)
%                                                 -(L*r*(5*Aya + 16*Ayb))/60
%                                                               -(Iyza*r)/10
%                                  -(r*(42*Izza + 7*Aa*L^2 + 15*Ab*L^2))/420];
%                  
% 
% M9 = [                                          (r*(3*Aya + 2*Ayb))/10
%                                                  (3*r*(Iyza + Iyzb))/(5*L)
%              (3*L*r*(3*Aa + 3*Ab))/140 - (3*r*(28*Iyya + 28*Iyyb))/(140*L)
%                                                   (L*r*(4*Aza + 5*Azb))/60
%                                   -(r*(6*Aa*L^2 - 42*Iyyb + 7*Ab*L^2))/420
%                                                                (Iyzb*r)/10
%                                                     (r*(2*Aya + 3*Ayb))/10
%                                                 -(3*r*(Iyza + Iyzb))/(5*L)
%                   (r*(21*Iyya + 21*Iyyb))/(35*L) + (L*r*(3*Aa + 10*Ab))/35
%                                                  (L*r*(5*Aza + 16*Azb))/60
%                                   (r*(42*Iyya + 7*Aa*L^2 + 15*Ab*L^2))/420
%                                                                (Iyza*r)/10];
%                                                
% M10 = [                           0
%               -(L*r*(5*Aya + 4*Ayb))/60
%                (L*r*(5*Aza + 4*Azb))/60
%                      (L*r*(Ja + Jb))/12
%                 -(L^2*r*(Aza + Azb))/60
%                 -(L^2*r*(Aya + Ayb))/60
%                                       0
%              -(L*r*(5*Aya + 16*Ayb))/60
%               (L*r*(5*Aza + 16*Azb))/60
%                    (L*r*(Ja + 3*Jb))/12
%                (L^2*r*(Aza + 2*Azb))/60
%                (L^2*r*(Aya + 2*Ayb))/60];
%            
%            
% M11 = [                    (L*r*(4*Aya + Ayb))/60
%                                           (Iyza*r)/10
%               (r*(7*Aa*L^2 - 42*Iyya + 6*Ab*L^2))/420
%                                (L^2*r*(Aza + Azb))/60
%  -(L*r*(14*Iyya + 14*Iyyb + 3*Aa*L^2 + 3*Ab*L^2))/840
%                               -(L*r*(Iyza + Iyzb))/60
%                                (L*r*(Aya - 6*Ayb))/60
%                                          -(Iyza*r)/10
%              (r*(42*Iyya + 7*Aa*L^2 + 15*Ab*L^2))/420
%                              (L^2*r*(Aza + 2*Azb))/60
%   (L*r*(28*Iyya + 84*Iyyb + 3*Aa*L^2 + 5*Ab*L^2))/840
%                              (L*r*(Iyza + 3*Iyzb))/30];
%                          
% 
% M12 = [                            (L*r*(4*Aza + Azb))/60
%                      -(r*(7*Aa*L^2 - 42*Izza + 6*Ab*L^2))/420
%                                                  -(Iyza*r)/10
%                                        (L^2*r*(Aya + Ayb))/60
%                                       -(L*r*(Iyza + Iyzb))/60
%          -(L*r*(14*Izza + 14*Izzb + 3*Aa*L^2 + 3*Ab*L^2))/840
%                                        (L*r*(Aza - 6*Azb))/60
%                     -(r*(42*Izza + 7*Aa*L^2 + 15*Ab*L^2))/420
%                                                   (Iyza*r)/10
%                                      (L^2*r*(Aya + 2*Ayb))/60
%                                      (L*r*(Iyza + 3*Iyzb))/30
%           (L*r*(28*Izza + 84*Izzb + 3*Aa*L^2 + 5*Ab*L^2))/840];
%       
% M = [M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12];

M = nan;

%% matrix transformation

kstr = Tr'*K*Tr; 
fstr = Tr'*F; 
mstr = Tr'*M*Tr; 

% Coord transformation matrix struct(x,y,z) to real(y,z,x)
struct2real = zeros(6,6);
struct2real(1,3) = 1;
struct2real(2,1) = 1;
struct2real(3,2) = 1;
struct2real(4,6) = 1;
struct2real(5,4) = 1;
struct2real(6,5) = 1;
s2r = zeros(12,12);
s2r(1:6, 1:6) = struct2real;
s2r(7:12, 7:12) = struct2real;
K2 = s2r * kstr * s2r';
F2 = s2r * fstr;
M2 = s2r * mstr * s2r';
