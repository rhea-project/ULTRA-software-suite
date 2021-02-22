function Wwing= Weight(AC,X,DV,DP)

ModeNum  = AC.Wing.Airfoils.ModeNum;
nT       = AC.Structure.nT;
yT       = AC.Structure.yT;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
yAirfoil = AC.Wing.Airfoils.yAirfoil;

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

T = X(1:4*nT);

CST      = AC.Wing.Airfoils.CST;
Cheby   = X(4*nT+1:4*nT+2*ModeNum*nAirfoil);
Cheby   = reshape(Cheby,nAirfoil,2*ModeNum);

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
 
Sw = sum((Geom(2:end,4)+Geom(1:end-1,4)).*(Geom(2:end,2)-Geom(1:end-1,2)));
    
Delta_x_Sweep = Geom(end,1)+Geom(end,4)*0.5-Geom(1,4)*0.5;    
Delta_y_Sweep = Geom(end,2);
Delta_Sweep = Delta_x_Sweep/Delta_y_Sweep;

Sweep_half = atan(Delta_Sweep);

yNode = linspace(0,1,nT);    

%%

for i = 1:nT
    tu = T(i); tl = T(nT+i); tfs = T(2*nT+i); trs = T(3*nT+i);

    eta_i = yT(i);

    Chord = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*(eta_i*b)/b1).*heaviside(heaviside(b1-(eta_i*b))-0.1)...
             + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*((eta_i*b)-b1)/b2).*heaviside(heaviside((eta_i*b)-b1)-0.6);

    CST_i = zeros(1,size(CST,2));
    for j=1:size(CST,1)-1
        CST_i = CST_i + (meshgrid(CST(j,:),1)+ (meshgrid(CST(j+1,:),1)-meshgrid(CST(j,:),1)).*...
            meshgrid((eta_i-yAirfoil(j))/((yAirfoil(j+1)-yAirfoil(j))),1:size(CST,2))').*meshgrid(heaviside(heaviside(eta_i-(yAirfoil(j)))-0.1),1:size(CST,2))'.*...
            meshgrid(heaviside(heaviside((yAirfoil(j+1))-eta_i)-0.6),1:size(CST,2))';
    end
    if eta_i == yAirfoil(end)
        CST_i(end,:) = CST(end,:);
    end


    Cheby_i = zeros(1,size(Cheby,2));
    for j=1:size(CST,1)-1
        Cheby_i = Cheby_i + (meshgrid(Cheby(j,:),1)+ (meshgrid(Cheby(j+1,:),1)-meshgrid(Cheby(j,:),1)).*...
            meshgrid((eta_i-yAirfoil(j))/((yAirfoil(j+1)-yAirfoil(j))),1:size(Cheby,2))').*meshgrid(heaviside(heaviside(eta_i-(yAirfoil(j)))-0.1),1:size(Cheby,2))'.*...
            meshgrid(heaviside(heaviside((yAirfoil(j+1))-eta_i)-0.6),1:size(Cheby,2))';
    end
    if eta_i == yAirfoil(end)
        Cheby_i(end,:) = Cheby(end,:);
    end
    

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

    
    x = linspace(0,1,50);
    x = fs+(rs-fs)*x;
    
    yu = CSTairfoil(CST_i(1:5),x); 
    yl = CSTairfoil(CST_i(6:10),x);        

    [~, yu] = Airfoil_Cheby(Cheby_i(1:length(Cheby_i)/2),x',yu',1);   
    [~, yl] = Airfoil_Cheby(Cheby_i(length(Cheby_i)/2+1:end),x',yl',2);


    x = x.*Chord;
    yu = yu.*Chord;
    yl = yl.*Chord;

    for j=1:length(x)-1
        dsu(j) = sqrt((x(j+1)-x(j))^2+(yu(j+1)-yu(j))^2);
        dsl(j) = sqrt((x(j+1)-x(j))^2+(yl(j+1)-yl(j))^2);
        dx(j) = (x(j)-x(j+1));

        xx(j)  = 0.5*(x(j+1)+x(j));
        yxu(j) = 0.5*(yu(j+1)+yu(j));
        yxl(j) = 0.5*(yl(j+1)+yl(j));
    end

    Su = sum(dsu);
    Sl = sum(dsl);
    Au(i) = Su*tu;
    Al(i) = Sl*tl;

    ylfs = yl(1);
    yufs = yu(1);
    hfs = yufs - ylfs - tu - tl;

    ylrs = yl(end);
    yurs = yu(end);
    hrs = yurs - ylrs - tu - tl;

    Afs(i) = hfs*tfs;
    Ars(i) = hrs*trs;

end

for i=1:nT-1
    Wu(i)  = AC.St.Mat_u.rho*0.5*(Au(i)+Au(i+1))*(yNode(i+1)-yNode(i))*b/cos(Sweep_half);
    Wl(i)  = AC.St.Mat_l.rho*0.5*(Al(i)+Al(i+1))*(yNode(i+1)-yNode(i))*b/cos(Sweep_half);
    Wfs(i) = AC.St.Mat_fs.rho*0.5*(Afs(i)+Afs(i+1))*(yNode(i+1)-yNode(i))*b/cos(Sweep_half);
    Wrs(i) = AC.St.Mat_rs.rho*0.5*(Ars(i)+Ars(i+1))*(yNode(i+1)-yNode(i))*b/cos(Sweep_half);
    Wt(i)  = Wu(i) + Wl(i) + Wfs(i) + Wrs(i);  

end

Wwb = 2*sum(Wt);

Wwing = 1.5*Wwb + 15*Sw;

end
