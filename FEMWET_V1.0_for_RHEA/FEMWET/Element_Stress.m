function [g_tension_u, g_compression_u, g_buckling_u, g_tension_l, g_compression_l, g_buckling_l, g_shear_fs, g_shear_rs, g_buckling_fs,g_buckling_rs, Au, Al, Afs, Ars, Xu, Yu, Xl, Yl, Xfs, Yfs, Xrs, Yrs,  ...
    dg_tension_du_u, dg_tension_du_l, dg_compression_du_u, dg_compression_du_l, dg_buckling_du_u, dg_buckling_du_l, dg_shear_du_fs, dg_shear_du_rs, dg_buckling_du_fs, dg_buckling_du_rs, ...
    dg_tension_dx_u, dg_tension_dx_l, dg_compression_dx_u, dg_compression_dx_l, dg_buckling_dx_u, dg_buckling_dx_l, dg_shear_dx_fs, dg_shear_dx_rs, dg_buckling_dx_fs, dg_buckling_dx_rs, dAu_dx, dAl_dx, dAfs_dx, dArs_dx] ...
    = Element_Stress(AC,U0,nNodea,Conc,ele,X,nT,yT,nAirfoil,yAirfoil,nSkin,nSpar,s2r,cNode0,xea,yea,dcNode,dXea,dYea,DV,DP)

global Adjoint SF
 
nNode = sum(nNodea);
ModeNum  = AC.Wing.Airfoils.ModeNum;


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


Tu = T(1:nT); Tl = T(nT+1:2*nT); Tfs = T(2*nT+1:3*nT); Trs = T(3*nT+1:4*nT);

% CSTu = T(4*nT+1:4*nT+5*nAirfoil);
% CSTl = T(4*nT+5*nAirfoil+1:4*nT+2*5*nAirfoil);
CST      = AC.Wing.Airfoils.CST;
Cheby   = T(4*nT+1:4*nT+2*ModeNum*nAirfoil);
Cheby   = reshape(Cheby,nAirfoil,2*ModeNum);

Cr     = T(4*nT+2*ModeNum*nAirfoil+1);
lambda = T(4*nT+2*ModeNum*nAirfoil+2);
b      = T(4*nT+2*ModeNum*nAirfoil+3);
Gamma  = T(4*nT+2*ModeNum*nAirfoil+4);

b1 = AC.Wing.kink*b;
b2 = b - b1;

%       x                      y          z                   chord                           twist                       
Geom = [0                      0          AC.Wing.Geom(1,3)   Cr                               ;                           
        b1*tan(Gamma)          b1         AC.Wing.Geom(2,3)   Cr-AC.Wing.kink*b*tan(Gamma)     ;
        b*tan(Gamma)           b          AC.Wing.Geom(3,3)   lambda*Cr                       ];

Sweep = atan((Geom(end,1)+ AC.Aero.Xs*Geom(end,4)-AC.Aero.Xs*Geom(1,4))/b);

cNode = [cNode0(1:nNode) cNode0(nNode+1:2*nNode) cNode0(2*nNode+1:3*nNode)];

%% element geometry

% [cNode, xea, yea] = Node_Generation(AC,T,nT,yT,nAirfoil,yAirfoil,nNode);

% yNode = linspace(0,1,nNode);
eta_ai = AC.Aileron.Geom(1); 
eta_af = AC.Aileron.Geom(2);
yN1 = linspace(0,eta_ai,nNodea(1));
yN2 = linspace(eta_ai,eta_af,nNodea(2)+1);
yN3 = linspace(eta_af,1,nNodea(3)+1);
yNode = [yN1 yN2(2:end) yN3(2:end)];

eta_1 = yNode(Conc(ele,1));  % start of the element
eta_2 = yNode(Conc(ele,2));  % end of the element
eta_i = 0.5*(eta_1+eta_2);

Chord = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*(eta_i*b)/b1).*heaviside(heaviside(b1-(eta_i*b))-0.1)...
         + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*((eta_i*b)-b1)/b2).*heaviside(heaviside((eta_i*b)-b1)-0.6);

Chord = Chord*cos(Sweep);     

for j=1:length(AC.St.Spars(:,1));
    if eta_i == AC.St.Spars(j,1)
        fs = AC.St.Spars(j,2);
        rs = AC.St.Spars(j,3);
        break;
    elseif eta_i > AC.St.Spars(j,1) && eta_i < AC.St.Spars(j+1,1)
        fs = AC.St.Spars(j,2) + (AC.St.Spars(j+1,2)-AC.St.Spars(j,2))*(eta_i-AC.St.Spars(j,1))/(AC.St.Spars(j+1,1)-AC.St.Spars(j,1));
        rs = AC.St.Spars(j,3) + (AC.St.Spars(j+1,3)-AC.St.Spars(j,3))*(eta_i-AC.St.Spars(j,1))/(AC.St.Spars(j+1,1)-AC.St.Spars(j,1));
        break;
    end
end
    
for j=1:nT
    if eta_i == yT(j)
        tu = Tu(j);
        tl = Tl(j);
        tfs = Tfs(j);
        trs = Trs(j);
        break;
    elseif eta_i > yT(j) && eta_i < yT(j+1)
        tu =  Tu(j) + (Tu(j+1)-Tu(j))*(eta_i-yT(j))/(yT(j+1)-yT(j));
        tl =  Tl(j) + (Tl(j+1)-Tl(j))*(eta_i-yT(j))/(yT(j+1)-yT(j));
        tfs =  Tfs(j) + (Tfs(j+1)-Tfs(j))*(eta_i-yT(j))/(yT(j+1)-yT(j));
        trs =  Trs(j) + (Trs(j+1)-Trs(j))*(eta_i-yT(j))/(yT(j+1)-yT(j));
        break;
    end
end


% for i=1:nAirfoil
%     if eta_i == yAirfoil(i)
%         CSTu = CSTu((i-1)*5+1:(i-1)*5+5);
%         CSTl = CSTl((i-1)*5+1:(i-1)*5+5);
%         break;
%     elseif eta_i > yAirfoil(i) && eta_i < yAirfoil(i+1)
%         CSTu = CSTu((i-1)*5+1:(i-1)*5+5) + (CSTu((i+1-1)*5+1:(i+1-1)*5+5)-CSTu((i-1)*5+1:(i-1)*5+5))*(eta_i-yAirfoil(i))/(yAirfoil(i+1)-yAirfoil(i));
%         CSTl = CSTl((i-1)*5+1:(i-1)*5+5) + (CSTl((i+1-1)*5+1:(i+1-1)*5+5)-CSTl((i-1)*5+1:(i-1)*5+5))*(eta_i-yAirfoil(i))/(yAirfoil(i+1)-yAirfoil(i));
%         break;
%     end
% end


CST_i = zeros(1,size(CST,2));
for i=1:size(CST,1)-1
    CST_i = CST_i + (meshgrid(CST(i,:),1)+ (meshgrid(CST(i+1,:),1)-meshgrid(CST(i,:),1)).*...
        meshgrid((eta_i-yAirfoil(i))/((yAirfoil(i+1)-yAirfoil(i))),1:size(CST,2))').*meshgrid(heaviside(heaviside(eta_i-(yAirfoil(i)))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1))-eta_i)-0.6),1:size(CST,2))';
end
if eta_i == yAirfoil(end)
    CST_i(end,:) = CST(end,:);
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



x = linspace(0,1,nSkin+1);
x = fs+(rs-fs)*x;

yu = CSTairfoil(CST_i(1:5),x); 
yl = CSTairfoil(CST_i(6:10),x);        

[~, yu] = Airfoil_Cheby(Cheby_i(1:length(Cheby)/2),x',yu',1);   
[~, yl] = Airfoil_Cheby(Cheby_i(length(Cheby)/2+1:end),x',yl',2);


x = x.*Chord;
yu = yu.*Chord;
yl = yl.*Chord;

Xea = xea*Chord;
Yea = yea*Chord;

%%

C1(1) = cNode(Conc(ele,1),2);   % coordinate of element start
C1(2) = cNode(Conc(ele,1),3);
C1(3) = cNode(Conc(ele,1),1);

C2(1) = cNode(Conc(ele,2),2);   % coordinate of element end
C2(2) = cNode(Conc(ele,2),3);
C2(3) = cNode(Conc(ele,2),1);

L = sqrt(sum((C2-C1).^2));

Coy   = [0 1 0];
Coz   = [0 0 1];

lox = (C2(1)-C1(1))/L;         % Find direction cosines and lambda matrix.
mox = (C2(2)-C1(2))/L;
nox = (C2(3)-C1(3))/L;

lam = [lox mox nox; Coy; Coz];
zerom = zeros(3,3);
Tr = [lam zerom zerom zerom          % Transformation matrix
      zerom lam zerom zerom
      zerom zerom lam zerom
      zerom zerom zerom lam];

U = Tr * U0;  % displacement vector in local coordinate

%% wingbox geometry

for i=1:length(x)-1
    dsu(i) = sqrt((x(i+1)-x(i))^2+(yu(i+1)-yu(i))^2);
    dsl(i) = sqrt((x(i+1)-x(i))^2+(yl(i+1)-yl(i))^2);
    dx(i) = (x(i)-x(i+1));
    
    xx(i)  = 0.5*(x(i+1)+x(i));
    yxu(i) = 0.5*(yu(i+1)+yu(i));
    yxl(i) = 0.5*(yl(i+1)+yl(i));
end

Su = sum(dsu);
Sl = sum(dsl);
Au = Su*tu;
Al = Sl*tl;

yufs = CSTairfoil(CST(1:5),[fs fs+0.01]); 
ylfs = CSTairfoil(CST(6:10),[fs fs+0.01]);        
[~, yufs] = Airfoil_Cheby(Cheby(1:length(Cheby)/2),[fs fs+0.01]',yufs',1);   
[~, ylfs] = Airfoil_Cheby(Cheby(length(Cheby)/2+1:end),[fs fs+0.01]',ylfs',2);
yufs = yufs(1)*Chord;
ylfs = ylfs(1)*Chord;

% ylfs = CSTairfoil(CSTl,fs)*Chord;
% yufs = CSTairfoil(CSTu,fs)*Chord;


hfs = yufs - ylfs - tu - tl;

% ylrs = CSTairfoil(CSTl,rs)*Chord;
% yurs = CSTairfoil(CSTu,rs)*Chord;

yurs = CSTairfoil(CST(1:5),[rs rs+0.01]); 
ylrs = CSTairfoil(CST(6:10),[rs rs+0.01]);        
[~, yurs] = Airfoil_Cheby(Cheby(1:length(Cheby)/2),[rs rs+0.01]',yurs',1);   
[~, ylrs] = Airfoil_Cheby(Cheby(length(Cheby)/2+1:end),[rs rs+0.01]',ylrs',2);
yurs = yurs(1)*Chord;
ylrs = ylrs(1)*Chord;


hrs = yurs - ylrs - tu - tl;

xfs = fs*Chord*ones(1,nSpar)+tfs/2;
xrs = rs*Chord*ones(1,nSpar)-trs/2;

nfs = 0:nSpar-1;
yxfs = ylfs+tl+(hfs)/(2*nSpar) + nfs*(hfs)/nSpar;
yxrs = ylrs+tl+(hrs)/(2*nSpar) + nfs*(hrs)/nSpar;

yun  = yxu - tu/2;
yln  = yxl + tl/2;

Afs = hfs*tfs;
Ars = hrs*trs;

%% coordinates of elements

Xu = xx' - Xea;
Xl = xx' - Xea;

Yu = yun' - Yea;
Yl = yln' - Yea;

Xfs = xfs' - Xea;
Xrs = xrs' - Xea;

Yfs = yxfs' - Yea;
Yrs = yxrs' - Yea;


%% Allowables


for i=1:length(AC.St.Mat_u.t)-1
    if tu == AC.St.Mat_u.t(i)
        Ftu_L_u  = AC.St.Mat_u.Ftu.L.A(i);
        Ftu_LT_u = AC.St.Mat_u.Ftu.LT.A(i);
        Fty_L_u  = AC.St.Mat_u.Fty.L.A(i);
        Fty_LT_u = AC.St.Mat_u.Fty.LT.A(i);
        Fcy_L_u  = AC.St.Mat_u.Fcy.L.A(i);
        Fcy_LT_u = AC.St.Mat_u.Fcy.LT.A(i);
        Fsu_u    = AC.St.Mat_u.Fsu(i);
        break;
    elseif tu > AC.St.Mat_u.t(i) && tu < AC.St.Mat_u.t(i+1)
        Ftu_L_u  = AC.St.Mat_u.Ftu.L.A(i) + (AC.St.Mat_u.Ftu.L.A(i+1)-AC.St.Mat_u.Ftu.L.A(i))*(tu-AC.St.Mat_u.t(i))/(AC.St.Mat_u.t(i+1)-AC.St.Mat_u.t(i));
        Ftu_LT_u = AC.St.Mat_u.Ftu.LT.A(i) + (AC.St.Mat_u.Ftu.LT.A(i+1)-AC.St.Mat_u.Ftu.LT.A(i))*(tu-AC.St.Mat_u.t(i))/(AC.St.Mat_u.t(i+1)-AC.St.Mat_u.t(i));
        Fty_L_u  = AC.St.Mat_u.Fty.L.A(i) + (AC.St.Mat_u.Fty.L.A(i+1)-AC.St.Mat_u.Fty.L.A(i))*(tu-AC.St.Mat_u.t(i))/(AC.St.Mat_u.t(i+1)-AC.St.Mat_u.t(i));
        Fty_LT_u = AC.St.Mat_u.Fty.LT.A(i) + (AC.St.Mat_u.Fty.LT.A(i+1)-AC.St.Mat_u.Fty.LT.A(i))*(tu-AC.St.Mat_u.t(i))/(AC.St.Mat_u.t(i+1)-AC.St.Mat_u.t(i));
        Fcy_L_u  = AC.St.Mat_u.Fcy.L.A(i) + (AC.St.Mat_u.Fcy.L.A(i+1)-AC.St.Mat_u.Fcy.L.A(i))*(tu-AC.St.Mat_u.t(i))/(AC.St.Mat_u.t(i+1)-AC.St.Mat_u.t(i));
        Fcy_LT_u = AC.St.Mat_u.Fcy.LT.A(i) + (AC.St.Mat_u.Fcy.LT.A(i+1)-AC.St.Mat_u.Fcy.LT.A(i))*(tu-AC.St.Mat_u.t(i))/(AC.St.Mat_u.t(i+1)-AC.St.Mat_u.t(i));
        Fsu_u    = AC.St.Mat_u.Fsu(i) + (AC.St.Mat_u.Fsu(i+1)-AC.St.Mat_u.Fsu(i))*(tu-AC.St.Mat_u.t(i))/(AC.St.Mat_u.t(i+1)-AC.St.Mat_u.t(i));
        break;
    end
end

for i=1:length(AC.St.Mat_l.t)-1
    if tl == AC.St.Mat_l.t(i)
        Ftu_L_l  = AC.St.Mat_l.Ftu.L.A(i);
        Ftu_LT_l = AC.St.Mat_l.Ftu.LT.A(i);
        Fty_L_l  = AC.St.Mat_l.Fty.L.A(i);
        Fty_LT_l = AC.St.Mat_l.Fty.LT.A(i);
        Fcy_L_l  = AC.St.Mat_l.Fcy.L.A(i);
        Fcy_LT_l = AC.St.Mat_l.Fcy.LT.A(i);
        Fsu_l    = AC.St.Mat_l.Fsu(i);
        break;
    elseif tl > AC.St.Mat_l.t(i) && tl < AC.St.Mat_l.t(i+1)
        Ftu_L_l  = AC.St.Mat_l.Ftu.L.A(i) + (AC.St.Mat_l.Ftu.L.A(i+1)-AC.St.Mat_l.Ftu.L.A(i))*(tl-AC.St.Mat_l.t(i))/(AC.St.Mat_l.t(i+1)-AC.St.Mat_l.t(i));
        Ftu_LT_l = AC.St.Mat_l.Ftu.LT.A(i) + (AC.St.Mat_l.Ftu.LT.A(i+1)-AC.St.Mat_l.Ftu.LT.A(i))*(tl-AC.St.Mat_l.t(i))/(AC.St.Mat_l.t(i+1)-AC.St.Mat_l.t(i));
        Fty_L_l  = AC.St.Mat_l.Fty.L.A(i) + (AC.St.Mat_l.Fty.L.A(i+1)-AC.St.Mat_l.Fty.L.A(i))*(tl-AC.St.Mat_l.t(i))/(AC.St.Mat_l.t(i+1)-AC.St.Mat_l.t(i));
        Fty_LT_l = AC.St.Mat_l.Fty.LT.A(i) + (AC.St.Mat_l.Fty.LT.A(i+1)-AC.St.Mat_l.Fty.LT.A(i))*(tl-AC.St.Mat_l.t(i))/(AC.St.Mat_l.t(i+1)-AC.St.Mat_l.t(i));
        Fcy_L_l  = AC.St.Mat_l.Fcy.L.A(i) + (AC.St.Mat_l.Fcy.L.A(i+1)-AC.St.Mat_l.Fcy.L.A(i))*(tl-AC.St.Mat_l.t(i))/(AC.St.Mat_l.t(i+1)-AC.St.Mat_l.t(i));
        Fcy_LT_l = AC.St.Mat_l.Fcy.LT.A(i) + (AC.St.Mat_l.Fcy.LT.A(i+1)-AC.St.Mat_l.Fcy.LT.A(i))*(tl-AC.St.Mat_l.t(i))/(AC.St.Mat_l.t(i+1)-AC.St.Mat_l.t(i));
        Fsu_l    = AC.St.Mat_l.Fsu(i) + (AC.St.Mat_l.Fsu(i+1)-AC.St.Mat_l.Fsu(i))*(tl-AC.St.Mat_l.t(i))/(AC.St.Mat_l.t(i+1)-AC.St.Mat_l.t(i));
        break;
    end
end

for i=1:length(AC.St.Mat_fs.t)-1
    if tfs == AC.St.Mat_fs.t(i)
        Ftu_L_fs  = AC.St.Mat_fs.Ftu.L.A(i);
        Ftu_LT_fs = AC.St.Mat_fs.Ftu.LT.A(i);
        Fty_L_fs  = AC.St.Mat_fs.Fty.L.A(i);
        Fty_LT_fs = AC.St.Mat_fs.Fty.LT.A(i);
        Fcy_L_fs  = AC.St.Mat_fs.Fcy.L.A(i);
        Fcy_LT_fs = AC.St.Mat_fs.Fcy.LT.A(i);
        Fsu_fs    = AC.St.Mat_fs.Fsu(i);
        break;
    elseif tfs > AC.St.Mat_fs.t(i) && tfs < AC.St.Mat_fs.t(i+1)
        Ftu_L_fs  = AC.St.Mat_fs.Ftu.L.A(i)  + (AC.St.Mat_fs.Ftu.L.A(i+1)-AC.St.Mat_fs.Ftu.L.A(i))*(tfs-AC.St.Mat_fs.t(i))/(AC.St.Mat_fs.t(i+1)-AC.St.Mat_fs.t(i));
        Ftu_LT_fs = AC.St.Mat_fs.Ftu.LT.A(i) + (AC.St.Mat_fs.Ftu.LT.A(i+1)-AC.St.Mat_fs.Ftu.LT.A(i))*(tfs-AC.St.Mat_fs.t(i))/(AC.St.Mat_fs.t(i+1)-AC.St.Mat_fs.t(i));
        Fty_L_fs  = AC.St.Mat_fs.Fty.L.A(i)  + (AC.St.Mat_fs.Fty.L.A(i+1)-AC.St.Mat_fs.Fty.L.A(i))*(tfs-AC.St.Mat_fs.t(i))/(AC.St.Mat_fs.t(i+1)-AC.St.Mat_fs.t(i));
        Fty_LT_fs = AC.St.Mat_fs.Fty.LT.A(i) + (AC.St.Mat_fs.Fty.LT.A(i+1)-AC.St.Mat_fs.Fty.LT.A(i))*(tfs-AC.St.Mat_fs.t(i))/(AC.St.Mat_fs.t(i+1)-AC.St.Mat_fs.t(i));
        Fcy_L_fs  = AC.St.Mat_fs.Fcy.L.A(i)  + (AC.St.Mat_fs.Fcy.L.A(i+1)-AC.St.Mat_fs.Fcy.L.A(i))*(tfs-AC.St.Mat_fs.t(i))/(AC.St.Mat_fs.t(i+1)-AC.St.Mat_fs.t(i));
        Fcy_LT_fs = AC.St.Mat_fs.Fcy.LT.A(i) + (AC.St.Mat_fs.Fcy.LT.A(i+1)-AC.St.Mat_fs.Fcy.LT.A(i))*(tfs-AC.St.Mat_fs.t(i))/(AC.St.Mat_fs.t(i+1)-AC.St.Mat_fs.t(i));
        Fsu_fs    = AC.St.Mat_fs.Fsu(i)    + (AC.St.Mat_fs.Fsu(i+1)-AC.St.Mat_fs.Fsu(i))*(tfs-AC.St.Mat_fs.t(i))/(AC.St.Mat_fs.t(i+1)-AC.St.Mat_fs.t(i));
        break;
    end
end

for i=1:length(AC.St.Mat_rs.t)-1
    if trs == AC.St.Mat_rs.t(i)
        Ftu_L_rs  = AC.St.Mat_rs.Ftu.L.A(i);
        Ftu_LT_rs = AC.St.Mat_rs.Ftu.LT.A(i);
        Fty_L_rs  = AC.St.Mat_rs.Fty.L.A(i);
        Fty_LT_rs = AC.St.Mat_rs.Fty.LT.A(i);
        Fcy_L_rs  = AC.St.Mat_rs.Fcy.L.A(i);
        Fcy_LT_rs = AC.St.Mat_rs.Fcy.LT.A(i);
        Fsu_rs    = AC.St.Mat_rs.Fsu(i);
        break;
    elseif trs > AC.St.Mat_rs.t(i) && trs < AC.St.Mat_rs.t(i+1)
        Ftu_L_rs  = AC.St.Mat_rs.Ftu.L.A(i)  + (AC.St.Mat_rs.Ftu.L.A(i+1)-AC.St.Mat_rs.Ftu.L.A(i))*(trs-AC.St.Mat_rs.t(i))/(AC.St.Mat_rs.t(i+1)-AC.St.Mat_rs.t(i));
        Ftu_LT_rs = AC.St.Mat_rs.Ftu.LT.A(i) + (AC.St.Mat_rs.Ftu.LT.A(i+1)-AC.St.Mat_rs.Ftu.LT.A(i))*(trs-AC.St.Mat_rs.t(i))/(AC.St.Mat_rs.t(i+1)-AC.St.Mat_rs.t(i));
        Fty_L_rs  = AC.St.Mat_rs.Fty.L.A(i)  + (AC.St.Mat_rs.Fty.L.A(i+1)-AC.St.Mat_rs.Fty.L.A(i))*(trs-AC.St.Mat_rs.t(i))/(AC.St.Mat_rs.t(i+1)-AC.St.Mat_rs.t(i));
        Fty_LT_rs = AC.St.Mat_rs.Fty.LT.A(i) + (AC.St.Mat_rs.Fty.LT.A(i+1)-AC.St.Mat_rs.Fty.LT.A(i))*(trs-AC.St.Mat_rs.t(i))/(AC.St.Mat_rs.t(i+1)-AC.St.Mat_rs.t(i));
        Fcy_L_rs  = AC.St.Mat_rs.Fcy.L.A(i)  + (AC.St.Mat_rs.Fcy.L.A(i+1)-AC.St.Mat_rs.Fcy.L.A(i))*(trs-AC.St.Mat_rs.t(i))/(AC.St.Mat_rs.t(i+1)-AC.St.Mat_rs.t(i));
        Fcy_LT_rs = AC.St.Mat_rs.Fcy.LT.A(i) + (AC.St.Mat_rs.Fcy.LT.A(i+1)-AC.St.Mat_rs.Fcy.LT.A(i))*(trs-AC.St.Mat_rs.t(i))/(AC.St.Mat_rs.t(i+1)-AC.St.Mat_rs.t(i));
        Fsu_rs    = AC.St.Mat_rs.Fsu(i)    + (AC.St.Mat_rs.Fsu(i+1)-AC.St.Mat_rs.Fsu(i))*(trs-AC.St.Mat_rs.t(i))/(AC.St.Mat_rs.t(i+1)-AC.St.Mat_rs.t(i));
        break;
    end
end

%% upper panel stress

% stresses
Eu = AC.St.Mat_u.E;
Gu = Eu/(2*(1+AC.St.Mat_u.nu));

sigma_u = Eu*((U(7)/L - U(1)/L) + Xu*(-6*U(3)/L^2 + 4*U(5)/L +6*U(9)/L^2 + 2*U(11)/L) + Yu*(6*U(2)/L^2+4*U(6)/L-6*U(8)/L^2+2*U(12)/L))*SF;
tau_u   = Gu*(-2*U(5)+Yu*(U(10)/L - U(4)/L))*SF;

% failure 
tension_u = sigma_u >= 0;
compression_u = sigma_u < 0;

tau_max_u = (Fty_L_u + Fty_LT_u + Fcy_L_u + Fcy_LT_u)/4 * 2*Fsu_u/(Ftu_L_u + Ftu_LT_u);

g_tension_u = zeros(size(sigma_u));
g_compression_u = zeros(size(sigma_u));
g_buckling_u = zeros(size(sigma_u));


if tension_u
    g_tension_u = sigma_u/Fty_L_u + (tau_u/tau_max_u).^2 - 1;
elseif compression_u
    F_buckling_u = AC.St.F*sqrt(-sigma_u*Au*Eu/(AC.St.RP*Chord*(rs-fs)));

    g_compression_u = -sigma_u/Fcy_L_u + (tau_u/tau_max_u).^2 - 1;
    g_buckling_u    =  -sigma_u./F_buckling_u - 1;
end


%% lower panel

% stresses
El = AC.St.Mat_l.E;
Gl = El/(2*(1+AC.St.Mat_l.nu));

sigma_l = El*((U(7)/L - U(1)/L) + Xl*(-6*U(3)/L^2 + 4*U(5)/L +6*U(9)/L^2 + 2*U(11)/L) + Yl*(6*U(2)/L^2+4*U(6)/L-6*U(8)/L^2+2*U(12)/L))*SF;
tau_l   = Gl*(-2*U(5) + Yl*(U(10)/L - U(4)/L))*SF;

% failure
tension_l = sigma_l >= 0;
compression_l = sigma_l < 0;

tau_max_l = (Fty_L_l + Fty_LT_l + Fcy_L_l + Fcy_LT_l)/4 * 2*Fsu_l/(Ftu_L_l + Ftu_LT_l);

g_tension_l = zeros(size(sigma_l));
g_compression_l = zeros(size(sigma_l));
g_buckling_l = zeros(size(sigma_l));

if tension_l
    g_tension_l = sigma_l/Fty_L_l + (tau_l/tau_max_l).^2 - 1;
elseif compression_l
    F_buckling_l = AC.St.F*sqrt(-sigma_l*Al*El/(AC.St.RP*Chord*(rs-fs)));
    
    g_compression_l = -sigma_l/Fcy_L_l + (tau_l/tau_max_l).^2 - 1;
    g_buckling_l = -sigma_l./F_buckling_l- 1;
end


%% front spar

% stresses
Efs = AC.St.Mat_fs.E;
Gfs = Efs/(2*(1+AC.St.Mat_fs.nu));

% sigma_fs = Efs*((U(7)/L - U(1)/L) + Xfs*(-6*U(3)/L^2 + 4*U(5)/L +6*U(9)/L^2 + 2*U(11)/L) + Yfs*(6*U(2)/L^2+4*U(6)/L-6*U(8)/L^2+2*U(12)/L))*SF;
tau_fs = Gfs*(-Xfs*(U(10)/L - U(4)/L))*SF;

% failure
    abf = AC.St.RP/hfs;
    if abf <=6 
        Ks_fs =  0.0631*abf^4 - 1.0154*abf^3 +  5.9681*abf^2 - 15.5251*abf +  21.1027;
    elseif abf > 6
        Ks_fs = 5.2;
    end

tau_buckling_fs = Ks_fs*Efs*(tfs/hfs)^2;
tau_max_fs = (Fty_L_fs + Fty_LT_fs + Fcy_L_fs + Fcy_LT_fs)/4 * 2*Fsu_fs/(Ftu_L_fs + Ftu_LT_fs);

g_shear_fs = abs(tau_fs)/tau_max_fs - 1;
g_buckling_fs = abs(tau_fs)/tau_buckling_fs - 1;


%% rear spar

% stresses
Ers = AC.St.Mat_rs.E;
Grs = Ers/(2*(1+AC.St.Mat_rs.nu));

% sigma_rs = Ers*((U(7)/L - U(1)/L) + Xrs*(-6*U(3)/L^2 + 4*U(5)/L +6*U(9)/L^2 + 2*U(11)/L) + Yrs*(6*U(2)/L^2+4*U(6)/L-6*U(8)/L^2+2*U(12)/L))*SF;
tau_rs = Grs*(-Xrs*(U(10)/L - U(4)/L))*SF;

% failure
    abr = AC.St.RP/hrs;
    if abr <=6 
        Ks_rs = 0.0631*abr^4 - 1.0154*abr^3 +  5.9681*abr^2 - 15.5251*abr +  21.1027;
    elseif abr >6
        Ks_rs = 5.2;
    end
    
tau_buckling_rs = Ks_rs*Ers*(trs/hrs)^2;
tau_max_rs = (Fty_L_rs + Fty_LT_rs + Fcy_L_rs + Fcy_LT_rs)/4 * 2*Fsu_rs/(Ftu_L_rs + Ftu_LT_rs);

g_shear_rs =    abs(tau_rs)/tau_max_rs - 1;
g_buckling_rs = abs(tau_rs)/tau_buckling_rs - 1;


%% Adjoint

if Adjoint ==1 && nargout > 14
    
    
    % ***** dsigma_du  & dtau_du ******
    
    dsigma_du_u = SF*Eu*[-1/L*ones(1,nSkin); Yu'*6/L^2; -Xu'*6/L^2; zeros(1,nSkin); Xu'*4/L; Yu'*4/L; 1/L*ones(1,nSkin); -Yu'*6/L^2; Xu'*6/L^2; zeros(1,nSkin); Xu'*2/L; Yu'*2/L];        
    dtau_du_u   = SF*Gu*[zeros(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); -Yu'/L; -2*ones(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); Yu'/L; zeros(1,nSkin); zeros(1,nSkin)];

    dsigma_du_l = SF*El*[-1/L*ones(1,nSkin); Yl'*6/L^2; -Xl'*6/L^2; zeros(1,nSkin); Xl'*4/L; Yl'*4/L; 1/L*ones(1,nSkin); -Yl'*6/L^2; Xl'*6/L^2; zeros(1,nSkin); Xl'*2/L; Yl'*2/L];    
    dtau_du_l   = SF*Gl*[zeros(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); -Yl'/L; -2*ones(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); zeros(1,nSkin); Yl'/L; zeros(1,nSkin); zeros(1,nSkin)];   
   
    dtau_du_fs  = SF*Gfs*[zeros(1,nSpar); zeros(1,nSpar); zeros(1,nSpar); Xfs'/L; zeros(1,nSpar);zeros(1,nSpar); zeros(1,nSpar); zeros(1,nSpar); zeros(1,nSpar); -Xfs'/L; zeros(1,nSpar); zeros(1,nSpar)];
    
    dtau_du_rs  = SF*Grs*[zeros(1,nSpar); zeros(1,nSpar); zeros(1,nSpar); Xrs'/L; zeros(1,nSpar);zeros(1,nSpar); zeros(1,nSpar); zeros(1,nSpar); zeros(1,nSpar); -Xrs'/L; zeros(1,nSpar); zeros(1,nSpar)];
    
    % ****** dg_du *****   
    
    dg_tension_du_u     = zeros(12,nSkin);
    dg_tension_du_l     = zeros(12,nSkin);
    dg_compression_du_u = zeros(12,nSkin);
    dg_compression_du_l = zeros(12,nSkin);
    dg_buckling_du_u    = zeros(12,nSkin);
    dg_buckling_du_l    = zeros(12,nSkin);
    

    dtau_dx_du_u = zeros(12,nSkin);
    dtau_dx_du_l = zeros(12,nSkin);
    
    for e = 1:nSkin
        dtau_dx_du_u(:,e) = 2/(tau_max_u^2)*tau_u(e)*dtau_du_u(:,e);
        dtau_dx_du_l(:,e) = 2/(tau_max_l^2)*tau_l(e)*dtau_du_l(:,e);
    end
        
    if tension_u
        dg_tension_du_u     = 1/Fty_L_u*dsigma_du_u + dtau_dx_du_u;      
    elseif compression_u
        dg_compression_du_u = -1/Fcy_L_u*dsigma_du_u + dtau_dx_du_u; 
        for e=1:nSkin
            dg_buckling_du_u(:,e) = -1/F_buckling_u(e)*dsigma_du_u(:,e) + ...
                +sigma_u(e)/(F_buckling_u(e)^2) *AC.St.F*sqrt(Au*Eu/(AC.St.RP*Chord*(rs-fs)))*(-1/(2*(-sigma_u(e))^(1/2)))*dsigma_du_u(:,e);
        end      
    end
        
    
    if tension_l
        dg_tension_du_l  = 1/Fty_L_l*dsigma_du_l + dtau_dx_du_l;
    elseif compression_l
        dg_compression_du_l  = 1/Fcy_L_l*-dsigma_du_l+ dtau_dx_du_l;
        for e=1:nSkin
             dg_buckling_du_l(:,e) = -1/F_buckling_l(e)*dsigma_du_l(:,e) + ...
                +sigma_l(e)/(F_buckling_l(e)^2) *AC.St.F*sqrt(Al*El/(AC.St.RP*Chord*(rs-fs)))*(-1/(2*(-sigma_l(e))^(1/2)))*dsigma_du_l(:,e);
        end
    end
 
    
    dg_shear_du_fs    = zeros(12,nSpar);
    dg_buckling_du_fs = zeros(12,nSpar);
    dg_shear_du_rs    = zeros(12,nSpar);
    dg_buckling_du_rs = zeros(12,nSpar);
    
    for e=1:nSpar
        dg_shear_du_fs(:,e)    = 1/tau_max_fs*dtau_du_fs(:,e)*tau_fs(e)/abs(tau_fs(e)); 
        dg_buckling_du_fs(:,e) = 1/tau_buckling_fs*dtau_du_fs(:,e)*tau_fs(e)/abs(tau_fs(e));

        dg_shear_du_rs(:,e)    = 1/tau_max_rs *dtau_du_rs(:,e)*tau_rs(e)/abs(tau_rs(e));
        dg_buckling_du_rs(:,e) = 1/tau_buckling_rs*dtau_du_rs(:,e)*tau_rs(e)/abs(tau_rs(e));
    end
    
    
        % **** coordinates change
    for e=1:nSkin
        dg_tension_du_u(:,e)     = s2r* (Tr'* dg_tension_du_u(:,e)); 
        dg_compression_du_u(:,e) = s2r* (Tr'* dg_compression_du_u(:,e)); 
        dg_buckling_du_u(:,e)    = s2r* (Tr'* dg_buckling_du_u(:,e));
        
        dg_tension_du_l(:,e)     = s2r* (Tr'* dg_tension_du_l(:,e)); 
        dg_compression_du_l(:,e) = s2r* (Tr'* dg_compression_du_l(:,e)); 
        dg_buckling_du_l(:,e)    = s2r* (Tr'* dg_buckling_du_l(:,e));   
    end
    
    for e = 1:nSpar
        dg_shear_du_fs(:,e)    = s2r* (Tr'*dg_shear_du_fs(:,e));
        dg_shear_du_rs(:,e)    = s2r* (Tr'*dg_shear_du_rs(:,e));
        dg_buckling_du_fs(:,e) = s2r* (Tr'*dg_buckling_du_fs(:,e));
        dg_buckling_du_rs(:,e) = s2r* (Tr'*dg_buckling_du_rs(:,e));
    end
    
    
    % ******** dg_dx ************
    [ad_g_tension_u, ad_g_compression_u, ad_g_buckling_u, ad_g_tension_l, ad_g_compression_l, ad_g_buckling_l, ad_g_shear_fs, ad_g_shear_rs, ad_g_buckling_fs,ad_g_buckling_rs, ad_Au, ad_Al, ad_Afs, ad_Ars] = ...
        Element_Stress(AC,U0,nNodea,Conc,ele,gradientinit(X),nT,yT,nAirfoil,yAirfoil,nSkin,nSpar,s2r,cNode0,xea,yea,dcNode,dXea,dYea,DV,DP);
    

    dg_tension_dx_u     = zeros(nSkin,length(X));
    dg_tension_dx_l     = zeros(nSkin,length(X));
    dg_compression_dx_u = zeros(nSkin,length(X));
    dg_compression_dx_l = zeros(nSkin,length(X));
    dg_buckling_dx_u    = zeros(nSkin,length(X));
    dg_buckling_dx_l    = zeros(nSkin,length(X));
    
    if tension_u
        dg_tension_dx_u = full(ad_g_tension_u.dx);
    elseif compression_u
        dg_compression_dx_u = full(ad_g_compression_u.dx);
        dg_buckling_dx_u = full(ad_g_buckling_u.dx);
    end
      
    if tension_l
        dg_tension_dx_l = full(ad_g_tension_l.dx);
    elseif compression_l
        dg_compression_dx_l = full(ad_g_compression_l.dx);
        dg_buckling_dx_l = full(ad_g_buckling_l.dx);
    end
    

    dg_shear_dx_fs = full(ad_g_shear_fs.dx);
    dg_shear_dx_rs = full(ad_g_shear_rs.dx);
    dg_buckling_dx_fs = full(ad_g_buckling_fs.dx);
    dg_buckling_dx_rs = full(ad_g_buckling_rs.dx);
    
    dAu_dx = full(ad_Au.dx);
    dAl_dx = full(ad_Al.dx);
    dAfs_dx = full(ad_Afs.dx);
    dArs_dx = full(ad_Ars.dx);
    
    % ******* d_dXea
    [ad_g_tension_u, ad_g_compression_u, ad_g_buckling_u, ad_g_tension_l, ad_g_compression_l, ad_g_buckling_l, ad_g_shear_fs, ad_g_shear_rs, ad_g_buckling_fs,ad_g_buckling_rs, ad_Au, ad_Al, ad_Afs, ad_Ars] = ...
        Element_Stress(AC,U0,nNodea,Conc,ele,X,nT,yT,nAirfoil,yAirfoil,nSkin,nSpar,s2r,cNode0,gradientinit(xea),yea,dcNode,dXea,dYea,DV,DP);
    
    if tension_u
        dg_tension_dx_u = dg_tension_dx_u + ad_g_tension_u.dx*dXea;
    elseif compression_u
        dg_compression_dx_u = dg_compression_dx_u  +ad_g_compression_u.dx*dXea;
        dg_buckling_dx_u = dg_buckling_dx_u + ad_g_buckling_u.dx*dXea;
    end
    
    if tension_l
        dg_tension_dx_l = dg_tension_dx_l + ad_g_tension_l.dx*dXea;
    elseif compression_l
        dg_compression_dx_l = dg_compression_dx_l + ad_g_compression_l.dx*dXea;
        dg_buckling_dx_l = dg_buckling_dx_l + ad_g_buckling_l.dx*dXea;
    end
    
    dg_shear_dx_fs = dg_shear_dx_fs + ad_g_shear_fs.dx*dXea;
    dg_shear_dx_rs = dg_shear_dx_rs + ad_g_shear_rs.dx*dXea;
    dg_buckling_dx_fs = dg_buckling_dx_fs + ad_g_buckling_fs.dx*dXea;
    dg_buckling_dx_rs = dg_buckling_dx_rs + ad_g_buckling_rs.dx*dXea;
    
%     dAu_dx = dAu_dx + ad_Au.dx*dXea;
%     dAl_dx = dAl_dx + ad_Al.dx*dXea;
%     dAfs_dx = dAfs_dx + ad_Afs.dx*dXea;
%     dArs_dx = dArs_dx + ad_Ars.dx*dXea;
    
    % *** d_dYea
     [ad_g_tension_u, ad_g_compression_u, ad_g_buckling_u, ad_g_tension_l, ad_g_compression_l, ad_g_buckling_l, ad_g_shear_fs, ad_g_shear_rs, ad_g_buckling_fs,ad_g_buckling_rs, ad_Au, ad_Al, ad_Afs, ad_Ars] = ...
        Element_Stress(AC,U0,nNodea,Conc,ele,X,nT,yT,nAirfoil,yAirfoil,nSkin,nSpar,s2r,cNode0,xea,gradientinit(yea),dcNode,dXea,dYea,DV,DP);
    
    if tension_u
        dg_tension_dx_u = dg_tension_dx_u + ad_g_tension_u.dx*dYea;
    elseif compression_u
        dg_compression_dx_u = dg_compression_dx_u  +ad_g_compression_u.dx*dYea;
        dg_buckling_dx_u = dg_buckling_dx_u + ad_g_buckling_u.dx*dYea;
    end
    
    if tension_l
        dg_tension_dx_l = dg_tension_dx_l + ad_g_tension_l.dx*dYea;
    elseif compression_l
        dg_compression_dx_l = dg_compression_dx_l + ad_g_compression_l.dx*dYea;
        dg_buckling_dx_l = dg_buckling_dx_l + ad_g_buckling_l.dx*dYea;
    end
    
%     dg_shear_dx_fs = dg_shear_dx_fs + ad_g_shear_fs.dx*dYea;
%     dg_shear_dx_rs = dg_shear_dx_rs + ad_g_shear_rs.dx*dYea;
%     dg_buckling_dx_fs = dg_buckling_dx_fs + ad_g_buckling_fs.dx*dYea;
%     dg_buckling_dx_rs = dg_buckling_dx_rs + ad_g_buckling_rs.dx*dYea;
    
%     dAu_dx = dAu_dx + ad_Au.dx*dYea;
%     dAl_dx = dAl_dx + ad_Al.dx*dYea;
%     dAfs_dx = dAfs_dx + ad_Afs.dx*dYea;
%     dArs_dx = dArs_dx + ad_Ars.dx*dYea;

    % *** d_dcNode
     [ad_g_tension_u, ad_g_compression_u, ad_g_buckling_u, ad_g_tension_l, ad_g_compression_l, ad_g_buckling_l, ad_g_shear_fs, ad_g_shear_rs, ad_g_buckling_fs,ad_g_buckling_rs, ad_Au, ad_Al, ad_Afs, ad_Ars] = ...
        Element_Stress(AC,U0,nNodea,Conc,ele,X,nT,yT,nAirfoil,yAirfoil,nSkin,nSpar,s2r,gradientinit(cNode0),xea,yea,dcNode,dXea,dYea,DV,DP);
    
    for c = 1:length(cNode0)
        
    if tension_u
        dg_tension_dx_u = dg_tension_dx_u + ad_g_tension_u.dx(:,c)*dcNode(c,:);
    elseif compression_u
        dg_compression_dx_u = dg_compression_dx_u  + ad_g_compression_u.dx(:,c)*dcNode(c,:);
        dg_buckling_dx_u = dg_buckling_dx_u + ad_g_buckling_u.dx(:,c)*dcNode(c,:);
    end
    
    if tension_l
        dg_tension_dx_l = dg_tension_dx_l + ad_g_tension_l.dx(:,c)*dcNode(c,:);
    elseif compression_l
        dg_compression_dx_l = dg_compression_dx_l + ad_g_compression_l.dx(:,c)*dcNode(c,:);
        dg_buckling_dx_l = dg_buckling_dx_l + ad_g_buckling_l.dx(:,c)*dcNode(c,:);
    end
    
    dg_shear_dx_fs = dg_shear_dx_fs + ad_g_shear_fs.dx(:,c)*dcNode(c,:);
    dg_shear_dx_rs = dg_shear_dx_rs + ad_g_shear_rs.dx(:,c)*dcNode(c,:);
    dg_buckling_dx_fs = dg_buckling_dx_fs + ad_g_buckling_fs.dx(:,c)*dcNode(c,:);
    dg_buckling_dx_rs = dg_buckling_dx_rs + ad_g_buckling_rs.dx(:,c)*dcNode(c,:);
    
%     dAu_dx = dAu_dx + ad_Au.dx(:,c)*dcNode(c,:);
%     dAl_dx = dAl_dx + ad_Al.dx(:,c)*dcNode(c,:);
%     dAfs_dx = dAfs_dx + ad_Afs.dx(:,c)*dcNode(c,:);
%     dArs_dx = dArs_dx + ad_Ars.dx(:,c)*dcNode(c,:);       


    end
    
end

end
