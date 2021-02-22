function Load = Loads(AC,AS,Wwing,Wfuel)


Poption = 0;


%%

q = AC.Aero.q;

y = linspace(0,1,50);


Af      = zeros(length(y),1);
L_aero  = zeros(length(y),1);
L_wing  = zeros(length(y),1);
L_pp    = zeros(length(y),1);
L_fuel  = zeros(length(y),1);


%% fuel tank volume
yf      = linspace(AC.St.FuelTank(1),AC.St.FuelTank(2),10);
Chordf  = interp1(AC.Wing.Geom(:,2),AC.Wing.Geom(:,4),yf*AC.Wing.Span/2);
fsf     = interp1(AC.St.Spars(:,1),AC.St.Spars(:,2),yf);
rsf     = interp1(AC.St.Spars(:,1),AC.St.Spars(:,3),yf);
A       = zeros(length(yf),1);

for i=1:length(yf)
    
    x = AC.Wing.Airfoil_coord(1,:,1).*Chordf(i);
    [yu, yl] = Airfoil_interp(AC.Wing.Airfoil_coord,AC.Wing.eta,yf(i));
    yu = yu.*Chordf(i);
    yl = yl.*Chordf(i);
        
    SPu = spline(x,yu);
    SPl = spline(x,yl);
        
    A(i) = quad(@(x)ppval(SPu,x),fsf(i)*Chordf(i),rsf(i)*Chordf(i)) - quad(@(x)ppval(SPl,x),fsf(i)*Chordf(i),rsf(i)*Chordf(i));
           
end

for i=1:length(y)
    if y(i)>= AC.St.FuelTank(1) && y(i)< AC.St.FuelTank(2)
        Af(i) = interp1(yf,A,y(i));
    else
        Af(i) = 0;
    end
end


SPa = spline(yf*AC.Wing.Span/2,A);
Vfuel = quad(@(yf)ppval(SPa,yf),yf(1)*AC.Wing.Span/2,yf(end)*AC.Wing.Span/2);
        
%% Loads

   
for i = 1:length(y)
    
    % aerodynamic
    L_aero(i) = interp1([0 AS.AVL.strip.Yst' AC.Wing.Span/2],[AS.AVL.strip.ccl(1) AS.AVL.strip.ccl' 0]*q...
        ,y(i)*AC.Wing.Span/2);
     
     
    % weight
    L_wing(i) = interp1(Wwing.y,-AC.Load.n_max*9.81*Wwing.w,y(i));
       
    % power plant
    if isfield(AC,'wingEngine')
         for j = 1:size(AC.WingEngine,1)
             if y(i)<=AC.WingEngine(j,1) && y(i+1)>=AC.WingEngine(j,1) 
                 L_pp(i)  = - AC.n_max*AC.WingEngine(j,2)*9.81;
             end
         end
    end

    % fuel
    if y(i)>= AC.St.FuelTank(1) && y(i)< AC.St.FuelTank(2)
           
        if y(i+1)>AC.St.FuelTank(2)
            y2 = AC.St.FuelTank(2);
            Af2 = A(end);
        else
            y2 = y(i+1);
            Af2 = Af(i+1);
        end
        
        L_fuel(i) = L_fuel(i+1) + ((Af2+Af(i))/2*(y2-y(i))*AC.Wing.Span/2)/Vfuel*-AC.Load.n_max *9.81*Wfuel;
        
    end
    
   
end
 

L = L_aero + L_pp + L_fuel; % + L_wing;    


T_aero = AS.AVL.strip.cm_c4*q*AC.Wing.MAC.*AS.AVL.strip.chord;

S_aero = zeros(length(y),1);
S_wing =  zeros(length(y),1);
S_fuel =  zeros(length(y),1);
S_pp =  zeros(length(y),1);


for i = length(y)-1:-1:1
    S_aero(i) = S_aero(i+1) + (L_aero(i)+L_aero(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
    S_wing(i) = S_wing(i+1) + (L_wing(i)+L_wing(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
    S_fuel(i) = S_fuel(i+1) + (L_fuel(i)+L_fuel(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
    S_pp(i)   = S_pp(i+1)   + (L_pp(i)+L_pp(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
end

S = S_aero + S_wing + S_fuel + S_pp;


M_aero = zeros(length(y),1);
M_wing =  zeros(length(y),1);
M_fuel =  zeros(length(y),1);
M_pp =  zeros(length(y),1);

for i = length(y)-1:-1:1
    M_aero(i) = M_aero(i+1) + (S_aero(i)+S_aero(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
    M_wing(i) = M_wing(i+1) + (S_wing(i)+S_wing(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
    M_fuel(i) = M_fuel(i+1) + (S_fuel(i)+S_fuel(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
    M_pp(i)   = M_pp(i+1)   + (S_pp(i)+S_pp(i+1))/2*(y(i+1)-y(i))*AC.Wing.Span/2;
end

M = M_aero + M_wing + M_fuel + M_pp;
 
T = interp1([0 AS.AVL.strip.Yst' AC.Wing.Span/2],[T_aero(1) T_aero' T_aero(end)],y*AC.Wing.Span/2);

if Poption == 1
    figure
    hold on
    plot(y,L_aero,'-b');
    plot(y,L_wing,'-r');
    plot(y,L_fuel,'-g');
    plot(y,L_pp,'-c');
    plot(y,L,'-k');
    legend('Aerodynamic loads','Wing weight loads','Fuel weight loads','Power plant weight loads','Total loads');
    hold off
end


Load.y = y;
Load.L.L = L;
Load.L.L_aero = L_aero;
Load.L.L_wing = L_wing;
Load.L.L_fuel = L_fuel;
Load.L.L_pp = L_pp;
Load.S = S;
Load.M = M;
Load.T = T;


