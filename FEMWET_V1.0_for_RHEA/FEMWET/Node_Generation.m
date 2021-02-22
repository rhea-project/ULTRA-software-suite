function [cNode, Xea, Zea] = Node_Generation(AC,X,DV,DP)

nT       = AC.Structure.nT;
ModeNum  = AC.Wing.Airfoils.ModeNum;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
nNode    = AC.FEM.nNodea;


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

Cr     = X(4*nT+2*ModeNum*nAirfoil+1);
lambda = X(4*nT+2*ModeNum*nAirfoil+2);
b      = X(4*nT+2*ModeNum*nAirfoil+3);
Gamma  = X(4*nT+2*ModeNum*nAirfoil+4);

b1 = AC.Wing.kink*b;
b2 = b - b1;

%       x                      y          z                   chord
Geom = [0                      0          AC.Wing.Geom(1,3)   Cr; 
        b1*tan(Gamma)          b1         AC.Wing.Geom(2,3)   Cr-AC.Wing.kink*b*tan(Gamma) ;
        b*tan(Gamma)           b          AC.Wing.Geom(3,3)   lambda*Cr];
   
eta_ai = AC.Aileron.Geom(1); 
eta_af = AC.Aileron.Geom(2);

        
[Xea, Zea] = Elastic_Axis(AC,X);


yN1 = linspace(0,eta_ai,nNode(1));
yN2 = linspace(eta_ai,eta_af,nNode(2)+1);
yN3 = linspace(eta_af,1,nNode(3)+1);

yNode = [yN1 yN2(2:end) yN3(2:end)];


yNode = yNode*b;
xLE = (Geom(1,1)+(Geom(2,1)-Geom(1,1))*yNode/b1).*heaviside(heaviside(b1-yNode)-0.1)...
    + (Geom(2,1)+(Geom(3,1)-Geom(2,1))*(yNode-b1)/b2).*heaviside(heaviside(yNode-b1)-0.6);
zLE = (Geom(1,3)+(Geom(2,3)-Geom(1,3))*yNode/b1).*heaviside(heaviside(b1-yNode)-0.1)...
    + (Geom(2,3)+(Geom(3,3)-Geom(2,3))*(yNode-b1)/b2).*heaviside(heaviside(yNode-b1)-0.6);
Chord = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*yNode/b1).*heaviside(heaviside(b1-yNode)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(yNode-b1)/b2).*heaviside(heaviside(yNode-b1)-0.6);

xNode = xLE + Chord*Xea;
zNode = zLE + Chord*Zea;

cNode = [xNode yNode zNode]';  % coordinates of nodes

