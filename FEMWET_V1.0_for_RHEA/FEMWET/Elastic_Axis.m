function [Xea, Yea] = Elastic_Axis(AC,T)

nT       = AC.Structure.nT;
yT       = AC.Structure.yT;
ModeNum  = AC.Wing.Airfoils.ModeNum;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
yAirfoil = AC.Wing.Airfoils.yAirfoil;
CST      = AC.Wing.Airfoils.CST;
% Cheby    = AC.Wing.Airfoils.Chebyshev;

Cheby   = T(4*nT+1:4*nT+2*ModeNum*nAirfoil);
Cheby   = reshape(Cheby,nAirfoil,2*ModeNum);

Cr      = T(4*nT+2*ModeNum*nAirfoil+1);
lambda = T(4*nT+2*ModeNum*nAirfoil+2);
b      = T(4*nT+2*ModeNum*nAirfoil+3);
Gamma  = T(4*nT+2*ModeNum*nAirfoil+4);

b1 = AC.Wing.kink*b;
b2 = b - b1;

%       x                      y          z                   chord
Geom = [0                      0          AC.Wing.Geom(1,3)   Cr; 
        b1*tan(Gamma)          b1         AC.Wing.Geom(2,3)   Cr-AC.Wing.kink*b*tan(Gamma) ;
        b*tan(Gamma)           b          AC.Wing.Geom(3,3)   lambda*Cr];


Sweep = atan((Geom(end,1)+ AC.Aero.Xs*Geom(end,4)-AC.Aero.Xs*Geom(1,4))/b);

Sx = 0;
Sy = 0;



%% airfoils
CSTi = zeros(length(yT),size(CST,2));
for i=1:size(CST,1)-1
    CSTi = CSTi + (meshgrid(CST(i,:),1:nT)+ (meshgrid(CST(i+1,:),1:nT)-meshgrid(CST(i,:),1:nT)).*...
        meshgrid((yT-yAirfoil(i))/((yAirfoil(i+1)-yAirfoil(i))),1:size(CST,2))').*meshgrid(heaviside(heaviside(yT-(yAirfoil(i)))-0.1),1:size(CST,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1))-yT)-0.6),1:size(CST,2))';
end
if yT(end) == yAirfoil(end)
    CSTi(end,:) = CST(end,:);
end

Chebyshevi = zeros(length(yT),size(Cheby,2));
for i=1:size(Cheby,1)-1
    Chebyshevi = Chebyshevi + (meshgrid(Cheby(i,:),1:nT)+ (meshgrid(Cheby(i+1,:),1:nT)-meshgrid(Cheby(i,:),1:nT)).*...
        meshgrid((yT-yAirfoil(i))/((yAirfoil(i+1)-yAirfoil(i))),1:size(Cheby,2))').*meshgrid(heaviside(heaviside(yT-(yAirfoil(i)))-0.1),1:size(Cheby,2))'.*...
        meshgrid(heaviside(heaviside((yAirfoil(i+1))-yT)-0.6),1:size(Cheby,2))';
end       
if yT(end) == yAirfoil(end)
    Chebyshevi(end,:) = Cheby(end,:);
end


Chord = (Geom(1,4)+(Geom(2,4)-Geom(1,4))*yT*b/b1).*heaviside(heaviside(b1-yT*b)-0.1)...
    + (Geom(2,4)+(Geom(3,4)-Geom(2,4))*(yT*b-b1)/b2).*heaviside(heaviside(yT*b-b1)-0.6);

Chord = Chord*cos(Sweep);

spar = zeros(length(yT),2);
for i=1:size(AC.St.Spars,1)-1
    spar = spar + (meshgrid(AC.St.Spars(i,2:3),1:nT)+ (meshgrid(AC.St.Spars(i+1,2:3),1:nT)-meshgrid(AC.St.Spars(i,2:3),1:nT)).*...
        meshgrid((yT-AC.St.Spars(i,1))/((AC.St.Spars(i+1,1)-AC.St.Spars(i,1))),1:2)').*meshgrid(heaviside(heaviside(yT-(AC.St.Spars(i,1)))-0.1),1:2)'.*...
        meshgrid(heaviside(heaviside((AC.St.Spars(i+1,1))-yT)-0.6),1:2)';
end       
if yT(end) == yAirfoil(end)
    spar(end,:) = AC.St.Spars(end,2:3);
end




for n=1:nT
        
   [xea, yea] =  Shear_center(CSTi(n,:),Chebyshevi(n,:),spar(n,1),spar(n,2),Chord(n),T(n),T(nT+n),T(2*nT+n),T(3*nT+n)); 
   
   Sx = Sx + xea;
   Sy = Sy + yea;
    
end

Xea = Sx/nT;
Yea = Sy/nT;

