function [Ltot, M, CL, Cl, cl2d, CCl, y, dy, CDi, dCDi_dWind] = Lift_VLM(X,AC,G,Wind,DV,DP)
%        [Ltot, M, Lift, y, dy, cl2d] = Lift_VLM(G,Mach,rho,Q,nNode,X,nT,nAirfoil,nChord,AC)

load Kcd.mat
Kcdi = Kcd(end,1);

G = G/sqrt(1-AC.Aero.Mach^2);  % Mach correction

nS       = sum(AC.Q3D.Grid.DnSpan);
nT       = AC.Structure.nT;
ModeNum  = AC.Wing.Airfoils.ModeNum;
nAirfoil = AC.Wing.Airfoils.nAirfoil;

%%


if DV == 1
    X  = [X(1:end-2); DP; X(end-1:end)];
elseif DV == 2
    X  = [DP(1:4*nT); X(1:end-2); DP(4*nT+1:end); X(end-1:end)];
elseif DV == 3 
    X  = [DP; X];
elseif DV ==4
    X  = [DP; X];
elseif DV ==5
    X  = [X(1:end-2); DP; X(end-1:end)];
elseif DV ==6
    X  = [X(1:4*nT); DP; X(4*nT+1:end)];
elseif DV ==7
    X  = X;
end

%%


Cr      = X(4*nT+2*ModeNum*nAirfoil+1);
lambda = X(4*nT+2*ModeNum*nAirfoil+2);
b      = X(4*nT+2*ModeNum*nAirfoil+3);
Gamma  = X(4*nT+2*ModeNum*nAirfoil+4);

b1 = AC.Wing.kink*b;
b2 = b - b1;

%       x                      y          z                   chord
Geom = [0                      0          AC.Wing.Geom(1,3)   Cr; 
        b1*tan(Gamma)          b1         AC.Wing.Geom(2,3)   Cr-AC.Wing.kink*b*tan(Gamma) ;
        b*tan(Gamma)           b          AC.Wing.Geom(3,3)   lambda*Cr];
           
Sweep = atan((Geom(end,1)+ AC.Aero.Xs*Geom(end,4)-AC.Aero.Xs*Geom(1,4))/b);


eta_ai = AC.Aileron.Geom(1); 
eta_af = AC.Aileron.Geom(2);

% Dy for lift calculation

yN1 = linspace(0,eta_ai,AC.Q3D.Grid.DnSpan(1));
yN2 = linspace(eta_ai,eta_af,AC.Q3D.Grid.DnSpan(2)+1);
yN3 = linspace(eta_af,1,AC.Q3D.Grid.DnSpan(3)+1);
yPanel = [yN1 yN2(2:end) yN3(2:end)]*b;

dy = yPanel(2:end) - yPanel(1:end-1);
Dy = dy';
for i=1:AC.Q3D.Grid.nChord-1
   Dy = [Dy;dy']; 
end


Gamma = G(nS:end)-G(1:end-nS+1);
Gamma = [G(1:nS-1); Gamma];

L = AC.Aero.rho*AC.Aero.Vinf*Gamma.*Dy;

%% TotalLift

Ltot = 2*sum(L);

Sw = sum((Geom(2:end,4)+Geom(1:end-1,4)).*(Geom(2:end,2)-Geom(1:end-1,2)));
CL = Ltot/(0.5*AC.Aero.rho*AC.Aero.Vinf^2*Sw);

%% Lift distribution

y = 0.5*(yPanel(2:end) + yPanel(1:end-1));
Lift = AC.Aero.rho*AC.Aero.Vinf*G(end-nS+2:end).*dy';
CCl = Lift./dy'/(0.5*AC.Aero.rho*AC.Aero.Vinf^2);

chord = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*y/b1).*heaviside(heaviside(b1-y)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(y-b1)/b2).*heaviside(heaviside(y-b1)-0.6);

Cl = CCl./chord';

%% Roling moment

M = sum(Lift.*y');

%% Induced drag

% % Wind = B*G;
% Wind = Wind/sqrt(1-AC.Aero.Mach^2);
% 
% Di = -AC.Aero.rho*Wind.*Gamma.*Dy;
% Ditot = 2*sum(Di);
% CDi = Ditot/(0.5*AC.Aero.rho*AC.Aero.Vinf^2*Sw);
% CDi = Kcdi(end)*CDi;
% 
% dCDi_dWind = Kcdi(end)*(1/(0.5*AC.Aero.rho*AC.Aero.Vinf^2*Sw)*2*ones(size(Di))).*(-AC.Aero.rho*Gamma.*Dy)/sqrt(1-AC.Aero.Mach^2);
% dCDi_dWind = dCDi_dWind'; 

% % Trefftz plane
Di = -AC.Aero.rho/2*(G(end-nS+2:end).*Wind.*dy'); 
CDi = Kcdi(end)*2*sum(Di)/(0.5*AC.Aero.rho*AC.Aero.Vinf^2*Sw);

dCDi_dWind =  Kcdi(end)*(1/(0.5*AC.Aero.rho*AC.Aero.Vinf^2*Sw)*2*ones(size(Di))).*(-AC.Aero.rho/2*G(end-nS+2:end).*dy');
dCDi_dWind = dCDi_dWind';

%%  for 2D sections
n2d      = AC.Q3D.N2D;
y2d      = linspace(0,1,n2d);
y2d      = y2d*(y(end)-y(1))+y(1);

cl2d = zeros(length(y2d),1);
for i=1:size(Cl,1)-1
    cl2d = cl2d + (meshgrid(Cl(i),1:n2d)+ (meshgrid(Cl(i+1),1:n2d)-meshgrid(Cl(i),1:n2d)).*...
        meshgrid((y2d-y(i))/((y(i+1)-y(i))),1:size(Cl,2))').*...
        meshgrid(heaviside(heaviside(y2d-(y(i)))-0.1),1:size(Cl,2))'.*...
        meshgrid(heaviside(heaviside((y(i+1))-y2d)-0.6),1:size(Cl,2))';
end
if y2d(end) == y(end)
    cl2d(end) = Cl(end);
end

cl2d = cl2d/cos(Sweep)^2;
