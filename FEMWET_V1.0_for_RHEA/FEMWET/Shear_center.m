function [xs, y0] = Shear_center(Coord,Cheby,fs,rs,Chord,tu,tl,tfs,trs)


CSTu = Coord(1:5);
CSTl = Coord(6:10);

x = linspace(fs,rs,50);

yu = CSTairfoil(CSTu,x); 
yl = CSTairfoil(CSTl,x);        

[~, yu] = Airfoil_Cheby(Cheby(1:length(Cheby)/2),x',yu',1);   
[~, yl] = Airfoil_Cheby(Cheby(length(Cheby)/2+1:end),x',yl',2);


x = x.*Chord;
yu = yu.*Chord;
yl = yl.*Chord;


%%
  
for i=1:length(x)-1
    dsu(i) = sqrt((x(i+1)-x(i))^2+(yu(i+1)-yu(i))^2);
    dsl(i) = sqrt((x(i+1)-x(i))^2+(yl(i+1)-yl(i))^2);
%     dx(i) = (x(i)-x(i+1));
    
    xx(i)  = 0.5*(x(i+1)+x(i));
    yxu(i) = 0.5*(yu(i+1)+yu(i));
    yxl(i) = 0.5*(yl(i+1)+yl(i));
    
    dydx_u(i) = (yu(i+1)-yu(i))/(x(i+1)-x(i));
    dydx_l(i) = (yl(i+1)-yl(i))/(x(i+1)-x(i));

end

yufs = yu(1);
ylfs = yl(1);
hfs = yufs-ylfs - tu - tl;

yurs = yu(end);
ylrs = yl(end);
hrs = yurs-ylrs - tu - tl;


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

y0   = (sum((yun*tu).*dsu)+sum((yln*tl).*dsl)+sum((yfs*tfs).*dsfs)+sum((yrs*trs).*dsrs))...
    /(sum(tu*dsu)+sum(tl*dsl)+sum(tfs*dsfs)+sum(trs*dsrs));

xs = 0.4*Chord;

yq = yq - y0;
xq0 = xq;
   

for iter = 1:100
    
    xq = xq0 - xs;

    tsq = [tu*ones(1,length(dsu)) tfs*ones(1,length(dsfs)) tl*ones(1,length(dsl)) trs*ones(1,length(dsrs))];

    Cu = yun - dydx_u.*xx;
    Cl = yln - dydx_l.*xxl;
    pu = abs(Cu./sqrt(1+dydx_u.^2));
    pl = abs(Cl./sqrt(1+dydx_l.^2));
    pfs = xfs;
    prs = xrs;

    P = [pu -pfs pl prs];

    Ixx = sum((dsq.*tsq).*yq.^2);
    Iyy = sum((dsq.*tsq).*xq.^2);
    Ixy = sum((dsq.*tsq).*(xq.*yq));

%     r = sqrt(xq.^2 + yq.^2);
%     J = sum((r.^2).*(dsq.*tsq));  % J as polar moment of inertia


    qbs = 0;
    qbs_plus = Ixy/(Ixx*Iyy-Ixy^2)*tsq.*xq.*dsq - Iyy/(Ixx*Iyy-Ixy^2)* tsq.*yq.*dsq;
    
    qs01 = 0;
    for s=1:length(dsq)
        qbs = qbs + qbs_plus(s);
        qs01 = qs01 + (qbs/tsq(s))*dsq(s);
    end
    
    qs0 = -qs01/sum(dsq./tsq);

    qbs = 0;
    xs_n = 0;
    for s=1:length(dsq)
        qbs = qbs + qbs_plus(s);
        qs = qbs + qs0;
        xs_n = xs_n + (qs*P(s))*dsq(s);
    end

    

    if abs(xs-xs_n)<1e-3
        break;
    else
        xs = xs_n;
    end
end

xs = xs/Chord;
y0 = y0/Chord;
