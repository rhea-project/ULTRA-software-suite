function [CL, CD, dCL_dX, dCD_dX, dCL_dAlpha, dCD_dAlpha, dCL_dAi, dCD_dAi, dCL_dU, dCD_dU, CDi, CDp, CDf, CDy, y] = Aero_coef(AC,X,Alpha,ai,U,iter,Par,cNode,dcNode_dX,DV,DP)


nX = length(X);
n2d = AC.Q3D.N2D;


%% step 1
XV = [X;Alpha;ai;U];
[da, dM, dR, dCST, dMode] = Section_aero(AC,gradientinit(XV),nX,DV,DP); 

aeff = da.x;   Meff = dM.x;   Reeff = dR.x;
da_dx = da.dx;  dM_dx = dM.dx;  dR_dx = dR.dx;
Mode = dMode.x;  dMode_dx = dMode.dx;  CSTi = dCST.x; % instead of dCST.dx

daeff_dX = da_dx(:,1:nX);    daeff_dAlpha = da_dx(:,nX+1);   daeff_dai = da_dx(:,nX+2:nX+1+n2d);  daeff_dU = da_dx(:,nX+2+n2d:end);
dMeff_dX = dM_dx(:,1:nX);    dMeff_dAlpha = dM_dx(:,nX+1);   dMeff_dai = dM_dx(:,nX+2:nX+1+n2d);  dMeff_dU = dM_dx(:,nX+2+n2d:end);
dReeff_dX = dR_dx(:,1:nX);   dReeff_dAlpha = dR_dx(:,nX+1);  dReeff_dai = dR_dx(:,nX+2:nX+1+n2d); dReeff_dU = dR_dx(:,nX+2+n2d:end);
dMode_dX = dMode_dx(:,1:nX);

%% step 2
Eff = [aeff; Meff; Reeff];

% [dCl, dCdp, dCdf] = Airfoil_analysis_emp(AC,n2d,gradientinit([Eff; Mode]));
% 
% Cleff = dCl.x;   Cdpeff = dCdp.x;      Cdfeff = dCdf.x;
% dCl_deff = dCl.dx;   dCdp_deff = dCdp.dx;    dCdf_deff = dCdf.dx;
% 
% dCleff_daeff = dCl_deff(:,1:n2d);  dCleff_dMeff = dCl_deff(:,n2d+1:2*n2d);   dCleff_dReeff = dCl_deff(:,2*n2d+1:3*n2d);
% dCdpeff_daeff = dCdp_deff(:,1:n2d);  dCdpeff_dMeff = dCdp_deff(:,n2d+1:2*n2d);   dCdpeff_dReeff = dCdp_deff(:,2*n2d+1:3*n2d);
% dCdfeff_daeff = dCdf_deff(:,1:n2d);  dCdfeff_dMeff = dCdf_deff(:,n2d+1:2*n2d);   dCdfeff_dReeff = dCdf_deff(:,2*n2d+1:3*n2d);
% dCleff_dMode = dCl_deff(:,3*n2d+1:end); dCdpeff_dMode = dCdp_deff(:,3*n2d+1:end); dCdfeff_dMode = dCdf_deff(:,3*n2d+1:end);

if Par ==0
    [Cleff, Cdpeff, Cdfeff, dCleff_daeff, dCleff_dMeff, dCleff_dReeff, dCdpeff_daeff, dCdpeff_dMeff, ...
    dCdpeff_dReeff, dCdfeff_daeff, dCdfeff_dMeff, dCdfeff_dReeff, dCleff_dMode, dCdpeff_dMode, dCdfeff_dMode] = ...
    Airfoil_analysis(AC,n2d,Eff,CSTi,Mode,iter);
elseif Par ==1
   [Cleff, Cdpeff, Cdfeff, dCleff_daeff, dCleff_dMeff, dCleff_dReeff, dCdpeff_daeff, dCdpeff_dMeff, ...
    dCdpeff_dReeff, dCdfeff_daeff, dCdfeff_dMeff, dCdfeff_dReeff, dCleff_dMode, dCdpeff_dMode, dCdfeff_dMode] = ...
    Airfoil_analysis_parallel(AC,n2d,Eff,CSTi,Mode,iter);
end
 
%% step 3
Cx = [X; Cleff; Cdpeff; Cdfeff; ai; U; cNode];

[dCL, dCD, dCDi, dCDp, dCDf,dCDy, dy] = Q3D_analysis(AC,gradientinit(Cx),nX,DV,DP);

CL = dCL.x;  CD = dCD.x; CDp = dCDp.x; CDf = dCDf.x; CDy = dCDy.x;  y = dy;  %y = dy.x;
CDi = dCDi; 
% CDi = dCDi.x;
dCL_dCx = dCL.dx;    dCD_dCx = dCD.dx;

dCL_dx = dCL_dCx(:,1:nX);  dCL_dCleff = dCL_dCx(:,nX+1:nX+n2d);  dCL_dCdpeff = dCL_dCx(:,nX+n2d+1:nX+2*n2d);
dCL_dCdfeff = dCL_dCx(:,nX+1+2*n2d:nX+3*n2d); dCL_dai = dCL_dCx(:,nX+3*n2d+1:nX+4*n2d); 
dCL_du = dCL_dCx(:,nX+4*n2d+1:nX+4*n2d+AC.FEM.nDOF);  dCL_dcNode = dCL_dCx(:,nX+4*n2d+AC.FEM.nDOF+1:end);
dCL_dx = dCL_dx + dCL_dcNode*dcNode_dX;

dCD_dx = dCD_dCx(:,1:nX);  dCD_dCleff = dCD_dCx(:,nX+1:nX+n2d);  dCD_dCdpeff = dCD_dCx(:,nX+n2d+1:nX+2*n2d); 
dCD_dCdfeff = dCD_dCx(:,nX+1+2*n2d:nX+3*n2d); dCD_dai = dCD_dCx(:,nX+3*n2d+1:nX+4*n2d);

% *******************
% dCD_dai = zeros(size(dCD_dai));
% *******************

%% step 4
for n=1:n2d
   % dCL_dx = []; if CL is not an explicit function of x

   dCL_dX(n,:) = dCL_dx(n,:) + dCL_dCleff(:,n)'*(dCleff_dMode*dMode_dX + dCleff_daeff*daeff_dX + dCleff_dMeff*dMeff_dX + dCleff_dReeff*dReeff_dX) ...
                         + dCL_dCdpeff(:,n)'*(dCdpeff_dMode*dMode_dX + dCdpeff_daeff*daeff_dX + dCdpeff_dMeff*dMeff_dX + dCdpeff_dReeff*dReeff_dX) ...
                         + dCL_dCdfeff(:,n)'*(dCdfeff_dMode*dMode_dX + dCdfeff_daeff*daeff_dX + dCdfeff_dMeff*dMeff_dX + dCdfeff_dReeff*dReeff_dX);


   dCL_dAlpha(n,:) =  dCL_dCleff(:,n)'*(dCleff_daeff*daeff_dAlpha + dCleff_dMeff*dMeff_dAlpha + dCleff_dReeff*dReeff_dAlpha) ...
                   + dCL_dCdpeff(:,n)'*(dCdpeff_daeff*daeff_dAlpha + dCdpeff_dMeff*dMeff_dAlpha + dCdpeff_dReeff*dReeff_dAlpha) ...
                   + dCL_dCdfeff(:,n)'*(dCdfeff_daeff*daeff_dAlpha + dCdfeff_dMeff*dMeff_dAlpha + dCdfeff_dReeff*dReeff_dAlpha);   


   dCL_dAi(n,:)  =  dCL_dai(n,:)  + dCL_dCleff(:,n)'*(dCleff_daeff*daeff_dai + dCleff_dMeff*dMeff_dai + dCleff_dReeff*dReeff_dai) ...
                       + dCL_dCdpeff(:,n)'*(dCdpeff_daeff*daeff_dai + dCdpeff_dMeff*dMeff_dai + dCdpeff_dReeff*dReeff_dai) ...
                       + dCL_dCdfeff(:,n)'*(dCdfeff_daeff*daeff_dai + dCdfeff_dMeff*dMeff_dai + dCdfeff_dReeff*dReeff_dai);  


   dCL_dU(n,:)  =  dCL_du(n,:) + dCL_dCleff(:,n)'*(dCleff_daeff*daeff_dU + dCleff_dMeff*dMeff_dU + dCleff_dReeff*dReeff_dU) ...
                               + dCL_dCdpeff(:,n)'*(dCdpeff_daeff*daeff_dU + dCdpeff_dMeff*dMeff_dU + dCdpeff_dReeff*dReeff_dU) ...
                               + dCL_dCdfeff(:,n)'*(dCdfeff_daeff*daeff_dU + dCdfeff_dMeff*dMeff_dU + dCdfeff_dReeff*dReeff_dU);  

end

dCD_dX = dCD_dx + dCD_dCleff*(dCleff_dMode*dMode_dX + dCleff_daeff*daeff_dX + dCleff_dMeff*dMeff_dX + dCleff_dReeff*dReeff_dX) ...
                   + dCD_dCdpeff*(dCdpeff_dMode*dMode_dX + dCdpeff_daeff*daeff_dX + dCdpeff_dMeff*dMeff_dX + dCdpeff_dReeff*dReeff_dX) ...
                   + dCD_dCdfeff*(dCdfeff_dMode*dMode_dX + dCdfeff_daeff*daeff_dX + dCdfeff_dMeff*dMeff_dX + dCdfeff_dReeff*dReeff_dX);

dCD_dAlpha =  dCD_dCleff*(dCleff_daeff*daeff_dAlpha + dCleff_dMeff*dMeff_dAlpha + dCleff_dReeff*dReeff_dAlpha) ...
                   + dCD_dCdpeff*(dCdpeff_daeff*daeff_dAlpha + dCdpeff_dMeff*dMeff_dAlpha + dCdpeff_dReeff*dReeff_dAlpha) ...
                   + dCD_dCdfeff*(dCdfeff_daeff*daeff_dAlpha + dCdfeff_dMeff*dMeff_dAlpha + dCdfeff_dReeff*dReeff_dAlpha);             

dCD_dAi  =  dCD_dai + dCD_dCleff*(dCleff_daeff*daeff_dai + dCleff_dMeff*dMeff_dai + dCleff_dReeff*dReeff_dai) ...
                       + dCD_dCdpeff*(dCdpeff_daeff*daeff_dai + dCdpeff_dMeff*dMeff_dai + dCdpeff_dReeff*dReeff_dai) ...
                       + dCD_dCdfeff*(dCdfeff_daeff*daeff_dai + dCdfeff_dMeff*dMeff_dai + dCdfeff_dReeff*dReeff_dai);     

% ********************
% dCD_dAi = zeros(size(dCD_dAi));
% *********************                   

dCD_dU =  dCD_dCleff*(dCleff_daeff*daeff_dU + dCleff_dMeff*dMeff_dU + dCleff_dReeff*dReeff_dU) ...
                       + dCD_dCdpeff*(dCdpeff_daeff*daeff_dU + dCdpeff_dMeff*dMeff_dU + dCdpeff_dReeff*dReeff_dU) ...
                       + dCD_dCdfeff*(dCdfeff_daeff*daeff_dU + dCdfeff_dMeff*dMeff_dU + dCdfeff_dReeff*dReeff_dU);   
                   
end
   
    



function [aeff, Meff, Reeff, CSTi, Chebyi] = Section_aero(AC,DV,nX,dv,DP)

n2d = AC.Q3D.N2D;
Vinf = AC.Aero.Vinf;
Mach = AC.Aero.Mach;

nT       = AC.Structure.nT;
ModeNum  = AC.Wing.Airfoils.ModeNum;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
yAirfoil = AC.Wing.Airfoils.yAirfoil;


X = DV(1:nX);
Alpha = DV(nX+1);
ai = DV(nX+2:nX+1+n2d);
U = DV(nX+2+n2d:end);


if dv == 1
    X  = [X(1:end-2); DP; X(end-1:end)];
elseif dv == 2
    X  = [DP(1:4*nT); X(1:end-2); DP(4*nT+1:end); X(end-1:end)];
elseif dv == 3 
    X  = [DP; X];
elseif dv ==4
    X  = [DP; X];
elseif dv ==5
    X  = [X(1:end-2); DP; X(end-1:end)];
elseif dv ==6
    X  = [X(1:4*nT); DP; X(4*nT+1:end)];
elseif dv ==7
    X  = X;
end

% Geometrical properties
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

% structural properties
eta_ai = AC.Aileron.Geom(1);
eta_af = AC.Aileron.Geom(2);
etaS1  = linspace(0,eta_ai,AC.Q3D.Grid.DnSpan(1));
etaS2  = linspace(eta_ai,eta_af,AC.Q3D.Grid.DnSpan(2)+1);
etaS3  = linspace(eta_af,1,AC.Q3D.Grid.DnSpan(3)+1);
yst    = [etaS1 etaS2(2:end) etaS3(2:end)]*b;
ys     = 0.5*(yst(2:end) + yst(1:end-1));
py     = U(5:6:end);

%% finding shape of 2D sections

y2d      = linspace(0,1,n2d);
% y2d      = y2d*(ys(end)-ys(1))+ys(1);
y2d      = y2d*b;

CST      = AC.Wing.Airfoils.CST;
Cheby   = X(4*nT+1:4*nT+2*ModeNum*nAirfoil);
Cheby   = reshape(Cheby,nAirfoil,2*ModeNum);


CSTi = zeros(length(y2d),size(CST,2));
for i=1:size(CST,1)-1
    CSTi = CSTi + (meshgrid(CST(i,:),1:n2d)+ (meshgrid(CST(i+1,:),1:n2d)-meshgrid(CST(i,:),1:n2d)).*...
        meshgrid((y2d-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(CST,2))').*...
        meshgrid(heaviside(heaviside(y2d-(yAirfoil(i)*b))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-y2d)-0.6),1:size(CST,2))';
end
if y2d(end) == yAirfoil(end)*b
    CSTi(end,:) = CST(end,:);
end


Chebyi = zeros(length(y2d),size(Cheby,2));
for i=1:size(Cheby,1)-1
    Chebyi = Chebyi + (meshgrid(Cheby(i,:),1:n2d)+ (meshgrid(Cheby(i+1,:),1:n2d)-meshgrid(Cheby(i,:),1:n2d)).*...
        meshgrid((y2d-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(Cheby,2))').*...
        meshgrid(heaviside(heaviside(y2d-(yAirfoil(i)*b))-0.1),1:size(Cheby,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-y2d)-0.6),1:size(Cheby,2))';
end
if y2d(end) == yAirfoil(end)*b
    Chebyi(end,:) = Cheby(end,:);
end
Chebyi = reshape(Chebyi',size(Chebyi,1)*size(Chebyi,2),1);



% airfoil chord
c2d = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*y2d/b1).*heaviside(heaviside(b1-y2d)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(y2d-b1)/b2).*heaviside(heaviside(y2d-b1)-0.6);

% structural twist
py2d = zeros(length(y2d),1);
for i=1:size(ys,2)-1
    py2d = py2d + (meshgrid(py(i),1:n2d)+ (meshgrid(py(i+1),1:n2d)-meshgrid(py(i),1:n2d)).*...
        meshgrid((y2d-yst(i))/((yst(i+1)-yst(i))),1:size(py,2))').*...
        meshgrid(heaviside(heaviside(y2d-(yst(i)))-0.1),1:size(py,2))'.*...
        meshgrid(heaviside(heaviside((yst(i+1))-y2d)-0.6),1:size(py,2))';
end
if y2d(end) == yst(end)
    py2d(end) = py(end);
end

% airfoil twist (aerodynamic+atructural
a2d = (Geom(1,5)+(Geom(2,5)-Geom(1,5))*y2d/b1).*heaviside(heaviside(b1-y2d)-0.1)...
    + (Geom(2,5)+(Geom(3,5)-Geom(2,5))*(y2d-b1)/b2).*heaviside(heaviside(y2d-b1)-0.6)  +  py2d';  % should be modified since py is arounf shear center


c2d    = c2d    /cos(Sweep);


%% Flow conditions

M2d  = Mach*cos(Sweep);          % section Mach
V2d  = Vinf*cos(Sweep);          % section speed

%% Effective 

aeff   = (Alpha + a2d')/cos(Sweep) - ai;
Meff   = M2d./cos(ai);
Veff   = V2d./cos(ai);
Reeff  = AC.Aero.rho*Veff.*c2d'/AC.Aero.mu; 

end

function [Cleff, Cdpeff, Cdfeff] = Airfoil_analysis_emp(AC,n2d,Eff)

ModeNum  = AC.Wing.Airfoils.ModeNum;

aeff = Eff(1:n2d);
Meff = Eff(n2d+1:2*n2d);
Reeff = Eff(2*n2d+1:3*n2d);
Mode = Eff(3*n2d+1:end);

Mode = reshape(Mode',2*ModeNum,n2d)';

%%

Cleff = 0.1*aeff*180/pi./sqrt(1-Meff.^2) + Mode(:,1) - 1.5*Mode(:,2);
Cdfeff = 0.455./(log(Reeff).^2.58.*(1+0.144*Meff.^2).^0.65);
Cdpeff = aeff;

end




function [CL, CD, CDi, CDp, CDf, CDy, y2d] = Q3D_analysis(AC,Cx,nX,DV,DP)

% global Kcdi Kcdp Kcdf
load Kcd
Kcdi = Kcd(end,1); Kcdp = Kcd(end,2); Kcdf = Kcd(end,3);  
% Kcdi = 1;   Kcdp = 1;   Kcdf = 1;

n2d = AC.Q3D.N2D;
nNode = AC.FEM.nNode;

X   = Cx(1:nX);
Cl  = Cx(nX+1:nX+n2d);
Cdp = Cx(nX+1+n2d:nX+2*n2d);
Cdf = Cx(nX+1+2*n2d:nX+3*n2d);
ai  = Cx(nX+1+3*n2d:nX+4*n2d);
U   = Cx(nX+1+4*n2d:nX+4*n2d+AC.FEM.nDOF);
cNode = Cx(nX+4*n2d+AC.FEM.nDOF+1:end);

px = U(4:6:end);
Ypx = cNode(nNode+1:2*nNode);


nT       = AC.Structure.nT;
ModeNum  = AC.Wing.Airfoils.ModeNum;
nAirfoil = AC.Wing.Airfoils.nAirfoil;


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

% Geometrical properties
Cr     = X(4*nT+2*ModeNum*nAirfoil+1);
lambda = X(4*nT+2*ModeNum*nAirfoil+2);
b      = X(4*nT+2*ModeNum*nAirfoil+3);
Gamma  = X(4*nT+2*ModeNum*nAirfoil+4);

b1 = AC.Wing.kink*b;
b2 = b - b1;

%       x                      y          z                   chord                           twist                       
Geom = [0                      0          AC.Wing.Geom(1,3)   Cr                               ;                           
        b1*tan(Gamma)          b1         AC.Wing.Geom(2,3)   Cr-AC.Wing.kink*b*tan(Gamma)     ;
        b*tan(Gamma)           b          AC.Wing.Geom(3,3)   lambda*Cr                       ];
 



Sweep = atan((Geom(end,1)+ AC.Aero.Xs*Geom(end,4)-AC.Aero.Xs*Geom(1,4))/b);
Sw = sum((Geom(2:end,4)+Geom(1:end-1,4)).*(Geom(2:end,2)-Geom(1:end-1,2)));

y2d      = linspace(0,1,n2d)*b;
dy2d = y2d(2:end) - y2d(1:end-1);

chord = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*y2d/b1).*heaviside(heaviside(b1-y2d)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(y2d-b1)/b2).*heaviside(heaviside(y2d-b1)-0.6);


px2d = zeros(length(y2d),size(px,2));
for i=1:size(px,1)-1
    px2d = px2d + (meshgrid(px(i),1:n2d)+ (meshgrid(px(i+1),1:n2d)-meshgrid(px(i),1:n2d)).*...
        meshgrid((y2d-Ypx(i))/((Ypx(i+1)-Ypx(i))),1:size(px,2))').*...
        meshgrid(heaviside(heaviside(y2d-(Ypx(i)))-0.1),1:size(px,2))'.*...
        meshgrid(heaviside(heaviside((Ypx(i+1))-y2d)-0.6),1:size(px,2))';
end
if y2d(end) == Ypx(end)
    px2d(end) = px(end);
end

%% 3D wing drag


cdp = Cdp./cos(ai) * cos(Sweep)^3;                     % pressure drag
cdf = Cdf./cos(ai);                                    % friction drag
cdi = zeros(size(cdf));  %Cl.*sin(ai)./cos(ai).^2 * cos(Sweep)^3;  
   

CDy = cdp+cdf+cdi;

CDi = sum(((cdi(1:end-1).*chord(1:end-1)') + (cdi(2:end).*chord(2:end)') )/2 .*dy2d');
CDp = sum(((cdp(1:end-1).*chord(1:end-1)') + (cdp(2:end).*chord(2:end)') )/2 .*dy2d');
CDf = sum(((cdf(1:end-1).*chord(1:end-1)') + (cdf(2:end).*chord(2:end)') )/2 .*dy2d');


CDi = Kcdi(end)*CDi*2/Sw;
CDf = Kcdf(end)*CDf*2/Sw;
CDp = Kcdp(end)*CDp*2/Sw;

CD  = CDi + CDp + CDf;


%% Section normal lift

CL = Cl.*cos(ai) - (Cdp+Cdf).*sin(ai);  % CL of the section perpendicular to the elastic axis

CL = CL.*cos(px2d);  % CL of the section perpendicular to the freestream

end