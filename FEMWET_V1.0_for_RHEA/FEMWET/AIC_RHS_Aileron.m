function [A, RHS, dRHS_dAlpha, B] = AIC_RHS_Aileron(X,AC,Alpha,U,Xea,DV,DP,delta_a)

nT       = AC.Structure.nT;
ModeNum  = AC.Wing.Airfoils.ModeNum;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
yAirfoil = AC.Wing.Airfoils.yAirfoil;
nNodea   = AC.FEM.nNodea;

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

Q = AC.Aero.Vinf;
   
 
%% geometry definition

nS = sum(AC.Q3D.Grid.DnSpan);

nC1 = AC.Q3D.Grid.DnChord(1);
nC2 = AC.Q3D.Grid.DnChord(2);
nC = nC1+nC2;

nP = (nS-1)*(nC-1);  % number of Panels


Cr     = X(4*nT+2*ModeNum*nAirfoil+1);
lambda = X(4*nT+2*ModeNum*nAirfoil+2);
b      = X(4*nT+2*ModeNum*nAirfoil+3);
Gamma  = X(4*nT+2*ModeNum*nAirfoil+4);
epsilon1  = X(4*nT+2*ModeNum*nAirfoil+5);
epsilon2  = X(4*nT+2*ModeNum*nAirfoil+6);

b1 = AC.Wing.kink*b;
b2 = b - b1;

%       x                      y          z                   chord                           twist                       
Geom = [0                      0          AC.Wing.Geom(1,3)   Cr                               0;                           
        b1*tan(Gamma)          b1         AC.Wing.Geom(2,3)   Cr-AC.Wing.kink*b*tan(Gamma)     epsilon1;
        b*tan(Gamma)           b          AC.Wing.Geom(3,3)   lambda*Cr                        epsilon2];
    
    
Sweep = atan((Geom(end,1)+ AC.Aero.Xs*Geom(end,4)-AC.Aero.Xs*Geom(1,4))/b);

% ** Aileron    
eta_ai = AC.Aileron.Geom(1);
eta_af = AC.Aileron.Geom(2);
cf_c   = AC.Aileron.Geom(3);
if ~exist('delta_a','var')
    delta_a = AC.Aileron.Geom(4)*pi/180;
end


CST      = AC.Wing.Airfoils.CST;
Cheby   = X(4*nT+1:4*nT+2*ModeNum*nAirfoil);
Cheby   = reshape(Cheby,nAirfoil,2*ModeNum);

% Deformations
us = U(1:6:end);    
vs = U(2:6:end);    
ws = U(3:6:end);
px = U(4:6:end);
py = U(5:6:end);


%% position of the Grid points

etaS1 = linspace(0,eta_ai,nNodea(1));
etaS2 = linspace(eta_ai,eta_af,nNodea(2)+1);
etaS3 = linspace(eta_af,1,nNodea(3)+1);
etaS = [etaS1 etaS2(2:end) etaS3(2:end)];

etaC1 = linspace(0,1-cf_c,nC1);
etaC2 = linspace(1-cf_c,1,nC2+1);
etaC = [etaC1 etaC2(2:end)];

yGrid = etaS*(b1+b2);
xLE = (Geom(1,1)+(Geom(2,1)-Geom(1,1))*yGrid/b1).*heaviside(heaviside(b1-yGrid)-0.1)...
    + (Geom(2,1)+(Geom(3,1)-Geom(2,1))*(yGrid-b1)/b2).*heaviside(heaviside(yGrid-b1)-0.6);
zLE = (Geom(1,3)+(Geom(2,3)-Geom(1,3))*yGrid/b1).*heaviside(heaviside(b1-yGrid)-0.1)...
    + (Geom(2,3)+(Geom(3,3)-Geom(2,3))*(yGrid-b1)/b2).*heaviside(heaviside(yGrid-b1)-0.6);


cGrid = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*yGrid/b1).*heaviside(heaviside(b1-yGrid)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(yGrid-b1)/b2).*heaviside(heaviside(yGrid-b1)-0.6);

aGrid = (Geom(1,5)+(Geom(2,5)-Geom(1,5))*yGrid/b1).*heaviside(heaviside(b1-yGrid)-0.1)...
    + (Geom(2,5)+(Geom(3,5)-Geom(2,5))*(yGrid-b1)/b2).*heaviside(heaviside(yGrid-b1)-0.6)  +  py';  % should be modified since py is arounf shear center

% airfoils
CSTi = zeros(length(yGrid),size(CST,2));
for i=1:size(CST,1)-1
    CSTi = CSTi + (meshgrid(CST(i,:),1:nS)+ (meshgrid(CST(i+1,:),1:nS)-meshgrid(CST(i,:),1:nS)).*...
        meshgrid((yGrid-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(CST,2))').*meshgrid(heaviside(heaviside(yGrid-(yAirfoil(i)*b))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-yGrid)-0.6),1:size(CST,2))';
end
if yGrid(end) == yAirfoil(end)*b
    CSTi(end,:) = CST(end,:);
end

Chebyshevi = zeros(length(yGrid),size(Cheby,2));
for i=1:size(Cheby,1)-1
    Chebyshevi = Chebyshevi + (meshgrid(Cheby(i,:),1:nS)+ (meshgrid(Cheby(i+1,:),1:nS)-meshgrid(Cheby(i,:),1:nS)).*...
        meshgrid((yGrid-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(Cheby,2))').*meshgrid(heaviside(heaviside(yGrid-(yAirfoil(i)*b))-0.1),1:size(Cheby,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-yGrid)-0.6),1:size(Cheby,2))';
end       
if yGrid(end) == yAirfoil(end)*b
    Chebyshevi(end,:) = Cheby(end,:);
end

xGrid = meshgrid(cGrid,1:nC).* meshgrid(etaC,1:nS)' + meshgrid(xLE,1:nC);
yGrid = meshgrid(yGrid,1:nC);
zGrid = meshgrid(zLE,1:nC);

dx = xGrid(2:end,:)-xGrid(1:end-1,:);
% dy = yGrid(:,2:end)-yGrid(:,1:end-1);
% S = (0.5*(dx(:,2:end)+dx(:,1:end-1))).*(0.5*(dy(2:end,:)+dy(1:end-1,:)));  % panel area

% ************* Aileron ****************

xHL = xLE + cGrid*(1-cf_c);

dXa = (xGrid - meshgrid(xHL,1:nC))* (1-cos(delta_a));
dXa = dXa.*heaviside(dXa).*heaviside(yGrid-eta_ai*(b1+b2)+1e-9).*heaviside(eta_af*(b1+b2)-yGrid+1e-9);
dZa = (xGrid - meshgrid(xHL,1:nC))* sin(delta_a);
dZa = dZa.*heaviside(dZa).*heaviside(yGrid-eta_ai*(b1+b2)+1e-9).*heaviside(eta_af*(b1+b2)-yGrid+1e-9);

%% position of the vortices @ 25% chord

xVrtx = xGrid + 0.25*[dx; dx(end,:)] + meshgrid(us,1:nC) - dXa;
yVrtx = yGrid  + meshgrid(vs,1:nC);

etaCamb = (xVrtx -  meshgrid(xLE,1:nC))./meshgrid(cGrid,1:nC);

yu = CSTairfoil(CSTi(:,1:5),etaCamb(1:end-1,1)); %.*meshgrid(cGrid,1:nC-1);
yl = CSTairfoil(CSTi(:,6:10),etaCamb(1:end-1,1)); %.*meshgrid(cGrid,1:nC-1);

[~, yu] = Airfoil_Cheby(Chebyshevi(:,1:length(Cheby(1,:))/2),etaCamb(1:end-1,1),yu,1);   
[~, yl] = Airfoil_Cheby(Chebyshevi(:,length(Cheby(1,:))/2+1:end),etaCamb(1:end-1,1),yl,2);  

% ****************
% the input airfoils are perpendicular to the sweep line so *cos(Sweep)
% the input airfoils are perpendicular to the elastic axis so /cos(px)
yu = yu*cos(Sweep)./meshgrid(cos(px),1:size(yu));
yl = yl*cos(Sweep)./meshgrid(cos(px),1:size(yl));   

yu = yu.*meshgrid(cGrid,1:nC-1);
yl = yl.*meshgrid(cGrid,1:nC-1);
camber = [(yu+yl)/2; zeros(1,nS)];

zVrtx = zGrid + camber  +  meshgrid(ws,1:nC) - dZa;  

[xVrtx, zVrtx] = twist(xVrtx,zVrtx,meshgrid(aGrid,1:nC),meshgrid(xGrid(1,:) + Xea*(xGrid(end,:)-xGrid(1,:)),1:nC),zVrtx);

% calculating the normal vector

xA = xVrtx(2:end,2:end) - xVrtx(1:end-1,1:end-1);
yA = yVrtx(2:end,2:end) - yVrtx(1:end-1,1:end-1);
zA = zVrtx(2:end,2:end) - zVrtx(1:end-1,1:end-1);

xB = xVrtx(1:end-1,2:end) - xVrtx(2:end,1:end-1);
yB = yVrtx(1:end-1,2:end) - yVrtx(2:end,1:end-1);
zB = zVrtx(1:end-1,2:end) - zVrtx(2:end,1:end-1);

ni = yA.*zB - yB.*zA;
nj = xB.*zA - xA.*zB;
nk = xA.*yB - xB.*yA;

Ln = sqrt(ni.^2 + nj.^2 + nk.^2);

ni = ni./Ln;
nj = nj./Ln;
nk = nk./Ln;

%% position of the collection points @ 75% chord

usc = 0.5*(us(1:end-1)+us(2:end));
vsc = 0.5*(vs(1:end-1)+vs(2:end));
wsc = 0.5*(ws(1:end-1)+ws(2:end));
pxc = 0.5*(px(1:end-1)+px(2:end));
pyc = 0.5*(py(1:end-1)+py(2:end));


etaSc = 0.5*(etaS(1:end-1)+etaS(2:end));
etaCc = etaC(1:end-1) + 0.75*(etaC(2:end)-etaC(1:end-1));

yColec = etaSc*(b1+b2);
xLEc = (Geom(1,1)+(Geom(2,1)-Geom(1,1))*yColec/b1).*heaviside(heaviside(b1-yColec)-0.1)...
    + (Geom(2,1)+(Geom(3,1)-Geom(2,1))*(yColec-b1)/b2).*heaviside(heaviside(yColec-b1)-0.6);
zLEc = (Geom(1,3)+(Geom(2,3)-Geom(1,3))*yColec/b1).*heaviside(heaviside(b1-yColec)-0.1)...
    + (Geom(2,3)+(Geom(3,3)-Geom(2,3))*(yColec-b1)/b2).*heaviside(heaviside(yColec-b1)-0.6);

cColec = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*yColec/b1).*heaviside(heaviside(b1-yColec)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(yColec-b1)/b2).*heaviside(heaviside(yColec-b1)-0.6);


aColec = (Geom(1,5)+(Geom(2,5)-Geom(1,5))*yColec/b1).*heaviside(heaviside(b1-yColec)-0.1)...
    + (Geom(2,5)+(Geom(3,5)-Geom(2,5))*(yColec-b1)/b2).*heaviside(heaviside(yColec-b1)-0.6)  +  pyc';


CSTi = zeros(length(yColec),size(CST,2));
for i=1:size(CST,1)-1
    CSTi = CSTi + (meshgrid(CST(i,:),1:nS-1)+ (meshgrid(CST(i+1,:),1:nS-1)-meshgrid(CST(i,:),1:nS-1)).*...
        meshgrid((yColec-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(CST,2))').*meshgrid(heaviside(heaviside(yColec-(yAirfoil(i)*b))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-yColec)-0.6),1:size(CST,2))';
end
if yColec(end) == yAirfoil(end)*b
    CSTi(end,:) = CST(end,:);
end

Chebyshevi = zeros(length(yColec),size(Cheby,2));
for i=1:size(Cheby,1)-1
    Chebyshevi = Chebyshevi + (meshgrid(Cheby(i,:),1:nS-1)+ (meshgrid(Cheby(i+1,:),1:nS-1)-meshgrid(Cheby(i,:),1:nS-1)).*...
        meshgrid((yColec-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(Cheby,2))').*meshgrid(heaviside(heaviside(yColec-(yAirfoil(i)*b))-0.1),1:size(Cheby,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-yColec)-0.6),1:size(Cheby,2))';
end       
if yColec(end) == yAirfoil(end)*b
    Chebyshevi(end,:) = Cheby(end,:);
end

yu = CSTairfoil(CSTi(:,1:5),etaCc); %.*meshgrid(cColec,etaCc);
yl = CSTairfoil(CSTi(:,6:10),etaCc); %.*meshgrid(cColec,etaCc);

[~, yu] = Airfoil_Cheby(Chebyshevi(:,1:length(Cheby(1,:))/2),etaCc,yu,1);  
[~, yl] = Airfoil_Cheby(Chebyshevi(:,length(Cheby(1,:))/2+1:end),etaCc,yl,2);   

yu = yu*cos(Sweep)./meshgrid(cos(pxc),1:size(yu));
yl = yl*cos(Sweep)./meshgrid(cos(pxc),1:size(yl));    % the input airfoils are perpendicular to the sweep line

yu = yu.*meshgrid(cColec,etaCc);
yl = yl.*meshgrid(cColec,etaCc);
camber = (yu+yl)/2;

xColec = meshgrid(cColec,1:nC-1).* meshgrid(etaCc,1:nS-1)' + meshgrid(xLEc,1:nC-1)  + meshgrid(usc,1:nC-1);
yColec = meshgrid(yColec,1:nC-1)   +  meshgrid(vsc,1:nC-1);
zColec = meshgrid(zLEc,1:nC-1) + camber   +  meshgrid(wsc,1:nC-1);


% ********** Aileron **************
xHLc = xLEc + cColec*(1-cf_c);
dXac = (xColec - meshgrid(xHLc,1:nC-1))* (1-cos(delta_a));
dXac = dXac.*heaviside(dXac).*heaviside(yColec-eta_ai*(b1+b2)-1e-9).*heaviside(eta_af*(b1+b2)-yColec+1e-9);
dZac = (xColec - meshgrid(xHLc,1:nC-1))* sin(delta_a);
dZac = dZac.*heaviside(dZac).*heaviside(yColec-eta_ai*(b1+b2)-1e-9).*heaviside(eta_af*(b1+b2)-yColec+1e-9);

xColec = xColec - dXac;
zColec = zColec - dZac;


[xColec, zColec] = twist(xColec,zColec,meshgrid(aColec,1:nC-1),meshgrid(xColec(1,:) + Xea*(xColec(end,:)-xColec(1,:)),1:nC-1),zColec);


%% Ring Vortices

Mxc = meshgrid(xColec');   Myc = meshgrid(yColec');    Mzc = meshgrid(zColec');
Mxc = Mxc';                Myc = Myc';                 Mzc = Mzc';


[~,Mxv1] = meshgrid(xColec,xVrtx(1:end-1,1:end-1)');  [~,Mxv2] = meshgrid(xColec,xVrtx(1:end-1,2:end)');
Mxv1 = Mxv1';                                         Mxv2 = Mxv2';
[~,Mxv3] = meshgrid(xColec,xVrtx(2:end,2:end)');      [~,Mxv4] = meshgrid(xColec,xVrtx(2:end,1:end-1)');
Mxv3 = Mxv3';                                         Mxv4 = Mxv4';

[~,Myv1] = meshgrid(yColec,yVrtx(1:end-1,1:end-1)');  [~,Myv2] = meshgrid(yColec,yVrtx(1:end-1,2:end)');
Myv1 = Myv1';                                         Myv2 = Myv2';
[~,Myv3] = meshgrid(yColec,yVrtx(2:end,2:end)');      [~,Myv4] = meshgrid(yColec,yVrtx(2:end,1:end-1)');
Myv3 = Myv3';                                          Myv4 = Myv4';

[~,Mzv1] = meshgrid(zColec,zVrtx(1:end-1,1:end-1)');  [~,Mzv2] = meshgrid(zColec,zVrtx(1:end-1,2:end)');
Mzv1 = Mzv1';                                         Mzv2 = Mzv2';
[~,Mzv3] = meshgrid(zColec,zVrtx(2:end,2:end)');      [~,Mzv4] = meshgrid(zColec,zVrtx(2:end,1:end-1)');
Mzv3 = Mzv3';                                          Mzv4 = Mzv4';

 
[u1r,v1r,w1r] = Vortex(Mxc,Myc,Mzc,Mxv1,Myv1,Mzv1,Mxv2,Myv2,Mzv2,1);
[u2r,v2r,w2r] = Vortex(Mxc,Myc,Mzc,Mxv2,Myv2,Mzv2,Mxv3,Myv3,Mzv3,1);
[u3r,v3r,w3r] = Vortex(Mxc,Myc,Mzc,Mxv3,Myv3,Mzv3,Mxv4,Myv4,Mzv4,1);
[u4r,v4r,w4r] = Vortex(Mxc,Myc,Mzc,Mxv4,Myv4,Mzv4,Mxv1,Myv1,Mzv1,1);

[u1l,v1l,w1l] = Vortex(Mxc,Myc,Mzc,Mxv1,-(Myv1),Mzv1,Mxv2,-(Myv2),Mzv2,-1);
[u2l,v2l,w2l] = Vortex(Mxc,Myc,Mzc,Mxv2,-(Myv2),Mzv2,Mxv3,-(Myv3),Mzv3,-1);
[u3l,v3l,w3l] = Vortex(Mxc,Myc,Mzc,Mxv3,-(Myv3),Mzv3,Mxv4,-(Myv4),Mzv4,-1);
[u4l,v4l,w4l] = Vortex(Mxc,Myc,Mzc,Mxv4,-(Myv4),Mzv4,Mxv1,-(Myv1),Mzv1,-1);

u1 = u1r+u1l;    v1 = v1r-v1l;      w1 = w1r+w1l;
u2 = u2r+u2l;    v2 = v2r-v2l;      w2 = w2r+w2l;
u3 = u3r+u3l;    v3 = v3r-v3l;      w3 = w3r+w3l;
u4 = u4r+u4l;    v4 = v4r-v4l;      w4 = w4r+w4l;

u = u1+u2+u3+u4;
v = v1+v2+v3+v4;
w = w1+w2+w3+w4;

% us = u2+u4;
% vs = v2+v4;
% ws = w2+w4;

%% ***************** wake ***********************

Mxcw = meshgrid(xColec',1:nS-1);    Mycw = meshgrid(yColec',1:nS-1);    Mzcw = meshgrid(zColec',1:nS-1);
Mxcw = Mxcw';                       Mycw = Mycw';                       Mzcw = Mzcw';

[~,Mxw1] = meshgrid(xColec,xVrtx(end,1:end-1)');  [~,Mxw2] = meshgrid(xColec,xVrtx(end,2:end)');
Mxw1 = Mxw1';                                      Mxw2 = Mxw2';
Mxw3  = inf;                                       Mxw4 = Mxw3;   
% Mxw3 = Mxw2 + ones(size(Mxw1))*50*(b1+b2);         Mxw4 = Mxw1 + ones(size(Mxw1))*50*(b1+b2);

[~,Myw1] = meshgrid(yColec,yVrtx(end,1:end-1)');  [~,Myw2] = meshgrid(xColec,yVrtx(end,2:end)');
Myw1 = Myw1';                                      Myw2 = Myw2';
Myw3 = Myw2;                                       Myw4 = Myw1;

[~,Mzw1] = meshgrid(zColec,zVrtx(end,1:end-1)');  [~,Mzw2] = meshgrid(zColec,zVrtx(end,2:end)');
Mzw1 = Mzw1';                                      Mzw2 = Mzw2';
Mzw3   = Mzw2;                                      Mzw4 = Mzw1;   
% Mzw3 = Mzw2 + ones(size(Mxw1))*50*(b1+b2)*tan(Alpha);   Mzw4 = Mzw1 + ones(size(Mxw1))*50*(b1+b2)*tan(Alpha); 
% Maw = meshgrid(zColec,aGrid(end,1:end-1)'); Maw = Maw';
% Mzw3 = Mzw2 + 50*(b1+b2)*tan(Alpha+Maw);           Mzw4 = Mzw1 +  50*(b1+b2)*tan(Alpha+Maw);  



[u1wr,v1wr,w1wr] = Vortex(Mxcw,Mycw,Mzcw,Mxw1,Myw1,Mzw1,Mxw2,Myw2,Mzw2,1);
[u2wr,v2wr,w2wr] = Vortex(Mxcw,Mycw,Mzcw,Mxw2,Myw2,Mzw2,Mxw3,Myw3,Mzw3,1);
[u3wr,v3wr,w3wr] = Vortex(Mxcw,Mycw,Mzcw,Mxw3,Myw3,Mzw3,Mxw4,Myw4,Mzw4,1);
[u4wr,v4wr,w4wr] = Vortex(Mxcw,Mycw,Mzcw,Mxw4,Myw4,Mzw4,Mxw1,Myw1,Mzw1,1);

[u1wl,v1wl,w1wl] = Vortex(Mxcw,Mycw,Mzcw,Mxw1,-(Myw1),Mzw1,Mxw2,-(Myw2),Mzw2,-1);
[u2wl,v2wl,w2wl] = Vortex(Mxcw,Mycw,Mzcw,Mxw2,-(Myw2),Mzw2,Mxw3,-(Myw3),Mzw3,-1);
[u3wl,v3wl,w3wl] = Vortex(Mxcw,Mycw,Mzcw,Mxw3,-(Myw3),Mzw3,Mxw4,-(Myw4),Mzw4,-1);
[u4wl,v4wl,w4wl] = Vortex(Mxcw,Mycw,Mzcw,Mxw4,-(Myw4),Mzw4,Mxw1,-(Myw1),Mzw1,-1);

u1w = u1wr+u1wl;    v1w = v1wr-v1wl;        w1w = w1wr+w1wl;
u2w = u2wr+u2wl;    v2w = v2wr-v2wl;        w2w = w2wr+w2wl;
u3w = u3wr+u3wl;    v3w = v3wr-v3wl;        w3w = w3wr+w3wl;
u4w = u4wr+u4wl;    v4w = v4wr-v4wl;        w4w = w4wr+w4wl;

uw = u1w+u2w+u3w+u4w;
vw = v1w+v2w+v3w+v4w;
ww = w1w+w2w+w3w+w4w;

% uws = u2w+u4w;
% vws = v2w+v4w;
% wws = w2w+w4w;

%%  **************** AIC ******************
u = u + [zeros(nP,nP-nS+1) uw];
v = v + [zeros(nP,nP-nS+1) vw];
w = w + [zeros(nP,nP-nS+1) ww];

% us = us + [zeros(nP,nP-nS+1) uws];
% vs = vs + [zeros(nP,nP-nS+1) vws];
% ws = ws + [zeros(nP,nP-nS+1) wws];

Mni = meshgrid(ni');   Mnj = meshgrid(nj');    Mnk = meshgrid(nk');
Mni = Mni';            Mnj = Mnj';             Mnk = Mnk';

A = u.*Mni + v.*Mnj + w.*Mnk;
% B = us.*Mni + vs.*Mnj + ws.*Mnk;

%% **************** Trefftz Plane *********************

yT = yColec(end,:);
zT = zColec(end,:);

yT1 = meshgrid(yT')';     yT2 = meshgrid(yT);  yT2 = -yT2.*(eye(size(yT2))-1);
zT1 = meshgrid(zT')';     zT2 = meshgrid(zT);

B1 = (yT1-yT2)./((zT1-zT2).^2 + (yT1-yT2).^2) - (yT1+yT2)./((zT1-zT2).^2 + (yT1+yT2).^2);
B2 = -((yT1+yT1)./((zT1-zT1).^2 + (yT1+yT1).^2));
B2 = B2.*eye(size(B2));

B = (B1+B2)*1/(2*pi);

%% **************** RHS ******************
Ni = reshape(ni',nP,1);
Nj = reshape(nj',nP,1);
Nk = reshape(nk',nP,1);

Qinf = [Q*cos(Alpha); 0; Q*sin(Alpha)];

RHS = Qinf(1)*Ni + Qinf(2)*Nj + Qinf(3)*Nk;
RHS = -RHS;

dRHS_dAlpha = sin(Alpha)*Q*Ni - cos(Alpha)*Q*Nk;


end

function y=CSTairfoil(A,x)


N1 = 0.5;
N2 = 1;

C = ((x.^N1)).*(1-x).^N2;

%% create Bernstein polynomial

% n = length(A(1,:));

% for v = 0:n-1
%     Sx(v+1,:) = nchoosek(n-1,v)*x.^v.*(1-x).^(n-1-v);
% end

Sx1 = nchoosek(4,0)*x.^0.*(1-x).^(4-0);
Sx2 = nchoosek(4,1)*x.^1.*(1-x).^(4-1);
Sx3 = nchoosek(4,2)*x.^2.*(1-x).^(4-2);
Sx4 = nchoosek(4,3)*x.^3.*(1-x).^(4-3);
Sx5 = nchoosek(4,4)*x.^4.*(1-x).^(4-4);

B1 = meshgrid(A(:,1),1:length(x))'.*meshgrid(Sx1,1:size(A,1));
B2 = meshgrid(A(:,2),1:length(x))'.*meshgrid(Sx2,1:size(A,1));
B3 = meshgrid(A(:,3),1:length(x))'.*meshgrid(Sx3,1:size(A,1));
B4 = meshgrid(A(:,4),1:length(x))'.*meshgrid(Sx4,1:size(A,1));
B5 = meshgrid(A(:,5),1:length(x))'.*meshgrid(Sx5,1:size(A,1));

yb = B1+B2+B3+B4+B5;

y = meshgrid(C,1:size(A,1)).*yb;
y = y';

end

function [X, Y] = twist(x,y,twist,x_twist,y_twist)

xt = x - x_twist;
yt = y - y_twist;

X = xt.*cos(-twist) - yt.*sin(-twist) + x_twist;
Y = xt.*sin(-twist) + yt.*cos(-twist) + y_twist;
end
