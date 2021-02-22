function Vf = fuel_volume(X,AC,DV,DP)

nS       = AC.Q3D.N2D;
nT       = AC.Structure.nT;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
yAirfoil = AC.Wing.Airfoils.yAirfoil;
ModeNum  = AC.Wing.Airfoils.ModeNum;

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


CST      = AC.Wing.Airfoils.CST;

Cheby   = X(4*nT+1:4*nT+2*ModeNum*nAirfoil);
Cheby   = reshape(Cheby,nAirfoil,2*ModeNum);

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
    

    
eta_i = linspace(0,1,nS);
eta_i = AC.St.FuelTank(1) + (AC.St.FuelTank(2)-AC.St.FuelTank(1))*eta_i;
eta_i = eta_i*b;

%% airfoil shape at each strip

CSTi = zeros(length(eta_i),size(CST,2));
for i=1:size(CST,1)-1
    CSTi = CSTi + (meshgrid(CST(i,:),1:nS)+ (meshgrid(CST(i+1,:),1:nS)-meshgrid(CST(i,:),1:nS)).*...
        meshgrid((eta_i-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(CST,2))').*meshgrid(heaviside(heaviside(eta_i-(yAirfoil(i)*b))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-eta_i)-0.6),1:size(CST,2))';
end
if eta_i(end) == yAirfoil(end)*b
    CSTi(end,:) = CST(end,:);
end

Chebyshevi = zeros(length(eta_i),size(Cheby,2));
for i=1:size(Cheby,1)-1
    Chebyshevi = Chebyshevi + (meshgrid(Cheby(i,:),1:nS)+ (meshgrid(Cheby(i+1,:),1:nS)-meshgrid(Cheby(i,:),1:nS)).*...
        meshgrid((eta_i-yAirfoil(i)*b)/((yAirfoil(i+1)-yAirfoil(i))*b),1:size(Cheby,2))').*meshgrid(heaviside(heaviside(eta_i-(yAirfoil(i)*b))-0.1),1:size(Cheby,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1)*b)-eta_i)-0.6),1:size(Cheby,2))';
end       
if eta_i(end) == yAirfoil(end)*b
    Chebyshevi(end,:) = Cheby(end,:);
end

Chord = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*eta_i/b1).*heaviside(heaviside(b1-eta_i)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(eta_i-b1)/b2).*heaviside(heaviside(eta_i-b1)-0.6);


%% spars position

for n=1:length(eta_i)
    for i=1:length(AC.St.Spars(:,1));
        if eta_i(n)/b == AC.St.Spars(i,1)
            fs_i(n) = AC.St.Spars(i,2);
            rs_i(n) = AC.St.Spars(i,3);
                    break;
        elseif eta_i(n)/b > AC.St.Spars(i,1) && eta_i(n)/b < AC.St.Spars(i+1,1)
            fs_i(n) = AC.St.Spars(i,2) + (AC.St.Spars(i+1,2)-AC.St.Spars(i,2))*(eta_i(n)/b-AC.St.Spars(i,1))/(AC.St.Spars(i+1,1)-AC.St.Spars(i,1));
            rs_i(n) = AC.St.Spars(i,3) + (AC.St.Spars(i+1,3)-AC.St.Spars(i,3))*(eta_i(n)/b-AC.St.Spars(i,1))/(AC.St.Spars(i+1,1)-AC.St.Spars(i,1));
                    break;
        end
    end
end


%%

x = linspace(0,1,50);

for i=1:nS
    xi = fs_i(i)+(rs_i(i)-fs_i(i))*x;
    
    yu = CSTairfoil(CSTi(i,1:5),xi); 
    yl = CSTairfoil(CSTi(i,6:10),xi);        

    [~, yu] = Airfoil_Cheby(Chebyshevi(i,1:length(Cheby)/2),x',yu,1);   
    [~, yl] = Airfoil_Cheby(Chebyshevi(i,length(Cheby)/2+1:end),x',yl,2);

    t = yu-yl;
    
    for j=1:length(xi)-1
        area(j) = (t(j) + t(j+1))/2*(xi(j+1)-xi(j))*Chord(i)^2;
    end
    Area(i) = sum(area);    
end

for i=1:nS-1
    Vfi(i) = (Area(i) + Area(i+1))/2*(eta_i(i+1)-eta_i(i));
end

Vf = sum(Vfi)*2;


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