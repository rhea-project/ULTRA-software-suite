function [Au, Al, Afs, Ars, A, Ax, Ay, Ixx, Iyy, Ixy, J, EA, EAx, EAy,  EIxx, EIyy, EIxy, GA, GAx, GAy, GJ, Area] ... 
    = Section_Properties(Coord,Cheby,fs,rs,Chord,Mat_u,Mat_l,Mat_fs,Mat_rs,tu,tl,tfs,trs,xs,y0)


CSTu = Coord(1:5);
CSTl = Coord(6:10);

% x = linspace(fs,rs,50);
x = linspace(0,1,50);
x = fs+(rs-fs)*x;

yu = CSTairfoil(CSTu,x); 
yl = CSTairfoil(CSTl,x);        

[~, yu] = Airfoil_Cheby(Cheby(1:length(Cheby)/2),x',yu',1);   
[~, yl] = Airfoil_Cheby(Cheby(length(Cheby)/2+1:end),x',yl',2);


x = x.*Chord;
yu = yu.*Chord;
yl = yl.*Chord;

t = yu-yl;

xs = xs*Chord;
y0 = y0*Chord;

%%
  
for i=1:length(x)-1
    dsu(i) = sqrt((x(i+1)-x(i))^2+(yu(i+1)-yu(i))^2);
    dsl(i) = sqrt((x(i+1)-x(i))^2+(yl(i+1)-yl(i))^2);
    dx(i) = (x(i)-x(i+1));
    
    xx(i)  = 0.5*(x(i+1)+x(i));
    yxu(i) = 0.5*(yu(i+1)+yu(i));
    yxl(i) = 0.5*(yl(i+1)+yl(i));
    
    dydx_u(i) = (yu(i+1)-yu(i))/(x(i+1)-x(i));
    dydx_l(i) = (yl(i+1)-yl(i))/(x(i+1)-x(i));
    
    area(i) = (t(i) + t(i+1))/2*(x(i+1)-x(i));
end

Area = sum(area);

yufs = yu(1);
ylfs = yl(1);
hfs = yufs-ylfs - tu - tl;

yurs = yu(end);
ylrs = yl(end);
hrs = yurs-ylrs - tu - tl;


Su = sum(dsu);
Sl = sum(dsl);

Au = Su*tu;
Al = Sl*tl;
Afs = hfs*tfs;
Ars = hrs*trs;

A = Au+Al+Afs+Ars;

%%  box area

A_box = 0;

for n=1:length(x)-1
    c1 = yu(n) - yl(n);
    c2 = yu(n+1) - yl(n+1);
    Ab = (c1+c2)/2*(x(n+1)-x(n));
    A_box = A_box + Ab;
end

%% spars discretization

xfs = fs*Chord*ones(1,20)+tfs/2;
xrs = rs*Chord*ones(1,20)-trs/2;

nfs = 0:19;
yfs = ylfs + tl + hfs/40 + nfs*hfs/20;
yrs = ylrs + tl + hrs/40 + nfs*hrs/20;

dsfs = hfs/20*ones(1,20);
dsrs = hrs/20*ones(1,20);

%% Shear and Bending moment centre

yun = yxu - tu/2;
yln = yxl + tl/2;

xxl = flipud(xx')';

dsq = [dsu dsfs dsl dsrs];
yq = [yun yfs yln yrs];
xq = [xx xfs xxl xrs];
tsq = [tu*ones(1,length(dsu)) tfs*ones(1,length(dsfs)) tl*ones(1,length(dsl)) trs*ones(1,length(dsrs))];
 
yq = yq - y0;
xq = xq - xs;

Ixx = sum((dsq.*tsq).*yq.^2);
Iyy = sum((dsq.*tsq).*xq.^2);
Ixy = sum((dsq.*tsq).*(xq.*yq));

r = sqrt(xq.^2 + yq.^2);
J = sum((r.^2).*(dsq.*tsq));  % J as polar moment of inertia

%% shear center

% y0   = (sum((yun*tu).*dsu)+sum((yln*tl).*dsl)+sum((yfs*tfs).*dsfs)+sum((yrs*trs).*dsrs))...
%     /(sum(tu*dsu)+sum(tl*dsl)+sum(tfs*dsfs)+sum(trs*dsrs));
% 
% xs = 0.4*Chord;
% 
% yq = yq - y0;
% xq0 = xq;
%    
% 
% for iter = 1:100
%     
%     xq = xq0 - xs;
% 
%     tsq = [tu*ones(1,length(dsu)) tfs*ones(1,length(dsfs)) tl*ones(1,length(dsl)) trs*ones(1,length(dsrs))];
% 
%     Cu = yun - dydx_u.*xx;
%     Cl = yln - dydx_l.*xxl;
%     pu = abs(Cu./sqrt(1+dydx_u.^2));
%     pl = abs(Cl./sqrt(1+dydx_l.^2));
%     pfs = xfs;
%     prs = xrs;
% 
%     P = [pu -pfs pl prs];
% 
%     Ixx = sum((dsq.*tsq).*yq.^2);
%     Iyy = sum((dsq.*tsq).*xq.^2);
%     Ixy = sum((dsq.*tsq).*(xq.*yq));
% 
%     r = sqrt(xq.^2 + yq.^2);
%     J = sum((r.^2).*(dsq.*tsq));  % J as polar moment of inertia
% 
% 
%     qbs = 0;
%     qbs_plus = Ixy/(Ixx*Iyy-Ixy^2)*tsq.*xq.*dsq - Iyy/(Ixx*Iyy-Ixy^2)* tsq.*yq.*dsq;
%     
%     qs01 = 0;
%     for s=1:length(dsq)
%         qbs = qbs + qbs_plus(s);
%         qs01 = qs01 + (qbs/tsq(s))*dsq(s);
%     end
%     
%     qs0 = -qs01/sum(dsq./tsq);
% 
%     qbs = 0;
%     xs_n = 0;
%     for s=1:length(dsq)
%         qbs = qbs + qbs_plus(s);
%         qs = qbs + qs0;
%         xs_n = xs_n + (qs*P(s))*dsq(s);
%     end
% 
%     
% 
%     if abs(xs-xs_n)<1e-3
%         break;
%     else
%         xs = xs_n;
%     end
% end


%% Stiffness calculation EA EAy EAz EIyy EIzz EIyz GA GAy GAz GJ

if nargout>2
    % EA, EI and GJ calculation
    dA = tsq.*dsq;
    Eq = [Mat_u.E*ones(1,length(dsu)) Mat_fs.E*ones(1,length(dsfs)) Mat_l.E*ones(1,length(dsl)) Mat_rs.E*ones(1,length(dsrs))];
    Gq = [Mat_u.E/(2*(1+Mat_u.nu))*ones(1,length(dsu)) Mat_fs.E/(2*(1+Mat_fs.nu))*ones(1,length(dsfs)) Mat_l.E/(2*(1+Mat_l.nu))*ones(1,length(dsl)) Mat_rs.E/(2*(1+Mat_rs.nu))*ones(1,length(dsrs))];

    % A = sum(dA);
    Ax = sum(yq.*dA);
    Ay = sum(xq.*dA);

    %     Stiff.y    = ys;
    EA   = sum(Eq.*dA);
    EAx  = sum(Eq.*yq.*dA);
    EAy  = sum(Eq.*xq.*dA);
    EIxx = sum(Eq.*dA.*yq.^2);
    EIyy = sum(Eq.*dA.*xq.^2);
    EIxy = sum(Eq.*dA.*xq.*yq);

    GA   = sum(Gq.*dA);   
    GAx  = sum(Gq.*yq.*dA);
    GAy  = sum(Gq.*xq.*dA);

    GJ = 4*A_box^2/sum(dsq./(tsq.*Gq));  % J based on thin wall hollow cylinder

end

%%

% figure
% hold on
% plot(x/Chord, yu/Chord);
% plot(x/Chord,yl/Chord);
% plot(x/Chord,(yu-tu)/Chord);
% plot(x/Chord,(yl+tl)/Chord);
% 
% plot(xfs/Chord,yfs/Chord);
% plot(xrs/Chord,yrs/Chord);
% plot((xfs-tfs/2)/Chord,yfs/Chord);
% plot((xrs+trs/2)/Chord,yrs/Chord);
% 
% plot(xs/Chord,y0/Chord,'or')
% axis equal;
