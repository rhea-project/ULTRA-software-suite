
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        %
%         Aircraft Lifting Surface       %
%            Weight Estimation           %
%              by: A. Elham              %
%           A.Elham@tudelft.nl           %
%                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% WWM.m is the run file of WWM tool
% to compile use: mcc -mgv WWM -a *.m

function EMWET(In)


if exist([In '.init'],'file')==0
    error([[In '.init'] 'dose not exist']);
end
if exist([In '.load'],'file')==0
    error([[In '.load'] ' dose not exist']);
end

[I AS] = ReadInput(In);

if I.Disp==1
clc

disp('*******************************************')
disp('EMWET Student Version 1.5')
disp('Copyright <C> 2010    Ali Elham, TU Delft')
disp('ATTENTION: Student vesion is a simplified version of EMWET.')
% disp('In order to have full version of AEWeight contact A.Elham@tudelft.nl');
disp('*******************************************')
end

Geo = Geometry(I);

[W S] = Weight_Estimation(I,Geo,AS);

% Write output weight and weight distribution file
fid = fopen([In '.weight'], 'wt');  

fprintf(fid,'Wing total weight(kg) %g\n',W.W_Wing);
fprintf(fid,'Wingbox weight(kg) %g\n',W.W_Wingbox);
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'   y/(b/2)   Chord[m]   tu[mm]     tl[mm]     tfs[mm]   trs[mm]  \n');

for i=1:length(S.y)
    SW(i,:) = [round(S.y(i)*100)/100,round(S.Chord(i)*100)/100,round(S.tu(i)*10000)/10,round(S.tl(i)*10000)/10,round(S.tfs(i)*10000)/10,round(S.trs(i)*10000)/10];
    fprintf(fid,'   %g       %g      %g        %g         %g        %g\n',SW(i,:));
end

fclose(fid);

quit
end

%%
function [I AS] = ReadInput(name)

fid = fopen([name '.init'], 'r');
init = textscan(fid, '%s');
fclose(fid);

Init = init{1};

I.Weight.MTOW =str2double(Init(1));
I.Weight.ZFW = str2double(Init(2));
I.n_max = str2double(Init(3));

I.Wing(1).Area =str2double(Init(4));
I.Wing(1).Span = str2double(Init(5));
I.Wing(1).SectionNumber = str2double(Init(6));

for i=1:str2double(Init(7))
    I.Wing(1).AirfoilPosition(i) = str2double(Init(6+2*i));
    I.Wing(1).AirfoilName(i) = Init(7+2*i);
    
    if exist([I.Wing(1).AirfoilName{i}],'file')==0
        error([[I.Wing(1).AirfoilName{i}] ' dose not exist']);
    end
    
end

j = 7+2*i+1;
    
for i=1:I.Wing(1).SectionNumber
    I.Wing(i).WingSection.Chord = str2double(Init(j+(i-1)*6));
    I.Wing(i).WingSection.Xle = str2double(Init(j+(i-1)*6+1));
    I.Wing(i).WingSection.Yle= str2double(Init(j+(i-1)*6+2));
    I.Wing(i).WingSection.Zle =str2double(Init(j+(i-1)*6+3));
    I.Wing(i).WingSection.FrontSparPosition =str2double(Init(j+(i-1)*6+4));
    I.Wing(i).WingSection.RearSparPosition = str2double(Init(j+(i-1)*6+5));
end

j = j+(i-1)*6+5+1;

I.WingFuelTank.Ystart = str2double(Init(j));
I.WingFuelTank.Yend = str2double(Init(j+1));

I.PP(1).WingEngineNumber   =  str2double(Init(j+2))*2;

if I.PP(1).WingEngineNumber~=0
    
    for i=1:I.PP(1).WingEngineNumber/2
        I.PP(i).EnginePosition = str2double(Init(j+2*i+1))*I.Wing(1).Span/2;
        I.PP(i).EngineWeight = str2double(Init(j+2*i+2));
    end

    j = j+2*i+2;
else
    j=j+2;
end

I.Material.Wing.UpperPanel.E = str2double(Init(j+1));
I.Material.Wing.UpperPanel.rho = str2double(Init(j+2));
I.Material.Wing.UpperPanel.Sigma_tensile = str2double(Init(j+3));
I.Material.Wing.UpperPanel.Sigma_compressive = str2double(Init(j+4));

I.Material.Wing.LowerPanel.E = str2double(Init(j+5));
I.Material.Wing.LowerPanel.rho = str2double(Init(j+6));
I.Material.Wing.LowerPanel.Sigma_tensile = str2double(Init(j+7));
I.Material.Wing.LowerPanel.Sigma_compressive = str2double(Init(j+8));

I.Material.Wing.FrontSpar.E = str2double(Init(j+9));
I.Material.Wing.FrontSpar.rho = str2double(Init(j+10));
I.Material.Wing.FrontSpar.Sigma_tensile = str2double(Init(j+11));
I.Material.Wing.FrontSpar.Sigma_compressive = str2double(Init(j+12));

I.Material.Wing.RearSpar.E = str2double(Init(j+13));
I.Material.Wing.RearSpar.rho = str2double(Init(j+14));
I.Material.Wing.RearSpar.Sigma_tensile = str2double(Init(j+15));
I.Material.Wing.RearSpar.Sigma_compressive = str2double(Init(j+16));


I.Structure.Wing.UpperPanelEfficiency=str2double(Init(j+17));
I.Structure.Wing.RibPitch =str2double(Init(j+18));

I.Disp = str2double(Init(j+19));

%%
fid = fopen([name '.load'], 'r');
load = textscan(fid, '%f %f %f');
fclose(fid);

AS.Y = load{1};
AS.L = load{2};
AS.M = load{3};

end
%%
function Geo = Geometry(I)

if I.Disp ==1
    disp( 'Geometry Generation ...')
end


%% Airfoil creation

% initializing
A0_w   = ones(2*5,1); 

for i=1:length(I.Wing(1).AirfoilName)
    [A.Wing(:,i) Wing.tc(i) Airfoil_coord.Wing(i,:,1) Airfoil_coord.Wing(i,:,2) Airfoil_coord.Wing(i,:,3)] =airfoilfit(I.Wing(1).AirfoilName{i},A0_w,1,0);
end



S_w    = zeros(I.Wing(1).SectionNumber-1,1);
S_mac  = zeros(I.Wing(1).SectionNumber-1,1);
S_Xmac = zeros(I.Wing(1).SectionNumber-1,1);
S_Ymac = zeros(I.Wing(1).SectionNumber-1,1);
S_Zmac = zeros(I.Wing(1).SectionNumber-1,1);
S_tc   = zeros(I.Wing(1).SectionNumber-1,1);
Geo.Wing.tc = zeros(I.Wing(1).SectionNumber,1);


for i = 1:I.Wing(1).SectionNumber
    Geo.Wing.tc(i)    = interp1(I.Wing(1).AirfoilPosition,Wing.tc,I.Wing(i).WingSection.Yle/I.Wing(I.Wing(1).SectionNumber).WingSection.Yle);
    Geo.Wing.Chord(i) = I.Wing(i).WingSection.Chord;
    Geo.Wing.Xle(i)   = I.Wing(i).WingSection.Xle;
    Geo.Wing.Yle(i)   = I.Wing(i).WingSection.Yle;
    Geo.Wing.Zle(i)   = I.Wing(i).WingSection.Zle;
    Geo.Wing.FS(i)    = I.Wing(i).WingSection.FrontSparPosition;
    Geo.Wing.RS(i)    = I.Wing(i).WingSection.RearSparPosition;
end


for i = 1:I.Wing(1).SectionNumber-1
    S_w(i) = (I.Wing(i).WingSection.Chord+I.Wing(i+1).WingSection.Chord)*...
        (I.Wing(i+1).WingSection.Yle-I.Wing(i).WingSection.Yle);
    S_mac(i) = (I.Wing(i).WingSection.Chord^2+I.Wing(i+1).WingSection.Chord^2)*...
        (I.Wing(i+1).WingSection.Yle-I.Wing(i).WingSection.Yle);
    S_Xmac(i) = (I.Wing(i).WingSection.Chord*I.Wing(i).WingSection.Xle+...
        I.Wing(i+1).WingSection.Chord*I.Wing(i+1).WingSection.Xle)*...
        (I.Wing(i+1).WingSection.Yle-I.Wing(i).WingSection.Yle);
    S_Ymac(i) = (I.Wing(i).WingSection.Chord*I.Wing(i).WingSection.Yle+...
        I.Wing(i+1).WingSection.Chord*I.Wing(i+1).WingSection.Yle)*...
        (I.Wing(i+1).WingSection.Yle-I.Wing(i).WingSection.Yle);
    S_Zmac(i) = (I.Wing(i).WingSection.Chord*I.Wing(i).WingSection.Zle+...
        I.Wing(i+1).WingSection.Chord*I.Wing(i+1).WingSection.Zle)*...
        (I.Wing(i+1).WingSection.Yle-I.Wing(i).WingSection.Yle);
    S_tc (i) = (Geo.Wing.tc(i)*I.Wing(i).WingSection.Chord + Geo.Wing.tc(i+1)*...
        I.Wing(i+1).WingSection.Chord)*(I.Wing(i+1).WingSection.Yle-I.Wing(i).WingSection.Yle);
        
    
end



if isfield(I.Wing(1),'Area')
    Geo.Wing.Area = I.Wing(1).Area;
else
    Geo.Wing.Area    =  sum(S_w);
end

Geo.Wing.mac     =  sum(S_mac)/Geo.Wing.Area;
Geo.Wing.Xmac    =  sum(S_Xmac)/Geo.Wing.Area;
Geo.Wing.Ymac    =  sum(S_Ymac)/Geo.Wing.Area;
Geo.Wing.Zmac    =  sum(S_Zmac)/Geo.Wing.Area;
Geo.Wing.tc_m    =  sum(S_tc)/Geo.Wing.Area;
Geo.Wing.Span    =  2*I.Wing(I.Wing(1).SectionNumber).WingSection.Yle;
Geo.Wing.AR      =  Geo.Wing.Span^2/Geo.Wing.Area;
Geo.Wing.Ct      =  I.Wing(I.Wing(1).SectionNumber).WingSection.Chord;
Geo.Wing.Cr      =  I.Wing(1).WingSection.Chord;
Geo.Wing.Taper   =  Geo.Wing.Ct/Geo.Wing.Cr;
Geo.Wing.Sweep_le=  atan((I.Wing(I.Wing(1).SectionNumber).WingSection.Xle -...
    I.Wing(1).WingSection.Xle)/(I.Wing(I.Wing(1).SectionNumber).WingSection.Yle-...
    I.Wing(1).WingSection.Yle));
Geo.Wing.Sweep   =  atan((I.Wing(I.Wing(1).SectionNumber).WingSection.Xle + 0.25*...
    I.Wing(I.Wing(1).SectionNumber).WingSection.Chord -I.Wing(1).WingSection.Xle-...
    0.25*I.Wing(1).WingSection.Chord)/(I.Wing(I.Wing(1).SectionNumber).WingSection.Yle-...
    I.Wing(1).WingSection.Yle));
Geo.Wing.Sweep_half=  atan((I.Wing(I.Wing(1).SectionNumber).WingSection.Xle + 0.5*...
    I.Wing(I.Wing(1).SectionNumber).WingSection.Chord -I.Wing(1).WingSection.Xle-...
    0.5*I.Wing(1).WingSection.Chord)/(I.Wing(I.Wing(1).SectionNumber).WingSection.Yle-...
    I.Wing(1).WingSection.Yle));
Geo.Wing.Airfoil_coord = Airfoil_coord.Wing;


if isfield(I.Structure.Wing,'RibPitch')
    Geo.Wing.RibPitch = I.Structure.Wing.RibPitch;
else
    Geo.Wing.RibPitch = 0.55*sqrt(Geo.Wing.Cr*Geo.Wing.tc(1));
end
    
end
%%


function [A tc x yu yl]=airfoilfit(Airfoil,A0,type,Poption)

global xe ye

warning off

Airfoil_adres=sprintf('%s%s%s',Airfoil);
airfoil=importdata(Airfoil_adres);

aua=A0(1:length(A0)/2);
ala=A0(length(A0)/2+1:length(A0));

X=airfoil(:,1);
Y=airfoil(:,2);

for j=1:length(X)-1
        if X(j+1)<X(j)
            x_u(j)=X(j);
            y_u(j)=Y(j);
        else
            x_l(j)=X(j);
            y_l(j)=Y(j);   
        end
end

xe=x_u;
ye=y_u;

OPToptions = optimset('Display','off');

if type==1
    Aua=fminunc(@Err_CST,aua,OPToptions);
    Eu=Err_CST(Aua);
elseif type==2
    Aua=fminunc(@Err_chebychev,aua,OPToptions);
    Eu=Err_chebychev(Aua);
end

xe=x_l;
ye=y_l;

if type==1
    Ala=fminunc(@Err_CST,ala,OPToptions);
    El=Err_CST(Ala);
elseif type==2
    Ala=fminunc(@Err_chebychev,ala,OPToptions);
    El=Err_chebychev(Ala);
end


Nx=100;
x=[1 1-sin(pi/2/Nx:pi/2/Nx:pi/2)];

if type ==1
    yu=CSTairfoil(Aua,x);
    yl=CSTairfoil(Ala,x);
elseif type==2
    yu=chebychev(Aua,x);
    yl=chebychev(Ala,x);
end

tc = max(yu-yl);

if Poption==1
    hold on
    plot(x,yu,'--r','LineWidth',2);
    plot(x_u,y_u,'-b');
    plot(x,yl,'--r','LineWidth',2);
    plot(x_l,y_l,'-b');
    xlabel('x/c');
    ylabel('y/c');
    if type==1
        legend('Fitted Airfoil -CST','Initial Airfoil','location','NorthEast');
    elseif type==2
        legend('Fitted Airfoil - Chebyshev','Initial Airfoil','location','NorthEast');
    end
    grid;
    hold off
end

A(1:length(A0)/2)=Aua;
A(length(A0)/2+1:length(A0))=Ala;

end


function rv=Err_CST(A)
global xe ye

n=length(xe);
e=0;
y=CSTairfoil(A,xe);

for i=1:n
    e=e+(ye(i)-y(i))^2;
end
rv=e;

end


function rv=Err_chebychev(A)
global xe ye

n=length(xe);
e=0;
y=chebychev(A,xe);

for i=1:n
    e=e+(ye(i)-y(i))^2;
end
rv=e;

end
%%
function y=CSTairfoil(A,x)


N1 = 0.5;
N2 = 1;

C = ((x.^N1)).*(1-x).^N2;

% create Bernstein polynomial

n = length(A);

for v = 0:n-1
    Sx(v+1,:) = nchoosek(n-1,v)*x.^v.*(1-x).^(n-1-v);
end


yb = zeros(1,length(x));

for i = 1:n
    yb(1,:) = yb(1,:) + A(i).*Sx(i,:);
end

y = C.*yb;

end
%%
function [yu yl x] = Airfoil_interp(Coord,Position,Y)

n = length(Position);
nx = length(Coord(1,:,1));

Xi = Coord(1,:,1);
Yi = Position;
Zui = zeros(nx,n);
Zli = zeros(nx,n);

for i = 1:n
    for j = 1:nx
        Zui(j,i) = Coord(i,j,2);
        Zli(j,i) = Coord(i,j,3);
    end
end

 x = Coord(1,:,1);
 yu = interp2(Yi,Xi,Zui,Y,x);
 yl = interp2(Yi,Xi,Zli,Y,x);
end
%%

function [W_wing St] = Weight_Estimation(I,Geo,AS)

if I.Disp ==1
disp('Weight Estimation  ...')
end
%% Geometry

y = linspace(0,1,round((Geo.Wing.Span/2)/Geo.Wing.RibPitch));
%y = AS.Y;
x = Geo.Wing.Airfoil_coord(1,:,1);
g = 9.81;

Chord = zeros(length(y)-1,1);fs = zeros(length(y)-1,1);rs = zeros(length(y)-1,1);
Au = zeros(length(y)-1,1);Al = zeros(length(y)-1,1);Afs = zeros(length(y)-1,1);
Ars = zeros(length(y)-1,1);Wu = zeros(length(y)-1,1);Wl = zeros(length(y)-1,1);
tsu = zeros(length(y)-1,1);tsl = zeros(length(y)-1,1);Wfs = zeros(length(y)-1,1);
Wrs = zeros(length(y)-1,1);tfs = zeros(length(y)-1,1);trs = zeros(length(y)-1,1);
W_box = zeros(length(y)-1,1);




%% WingBox weight

Wwing.y = [0 1];
W_w0 = 0;%0.12*I.Weight.MTOW;
Wfuel = (I.Weight.MTOW - I.Weight.ZFW)*0.5;

for iter = 1:20
    
    if I.Disp==1   
        disp(' ')
        disp([' iteration      ' num2str(iter) '...'])
    end
    
    Wwing.w = [2*W_w0/Geo.Wing.Span 0];
    Load = Loads(I,Geo,AS,Wwing,Wfuel);


   for i = 1 :length(y)-1
        
        [yu yl] = Airfoil_interp(Geo.Wing.Airfoil_coord,I.Wing(1).AirfoilPosition,y(i));

        Coord(:,1) = x;
        Coord(:,2) = yu/cos(Geo.Wing.Sweep_half);
        Coord(:,3) = yl/cos(Geo.Wing.Sweep_half);

        Chord(i) = interp1(Geo.Wing.Yle/(Geo.Wing.Span/2),Geo.Wing.Chord,y(i))*cos(Geo.Wing.Sweep_half);
        fs(i)    = interp1(Geo.Wing.Yle/(Geo.Wing.Span/2),Geo.Wing.FS,y(i));
        rs(i)    = interp1(Geo.Wing.Yle/(Geo.Wing.Span/2),Geo.Wing.RS,y(i));

        [Au(i) Al(i) Afs(i) Ars(i) tsu(i) tsl(i) tfs(i) trs(i) Su(i) Sl(i) h_fs(i) h_rs(i)] = Section_Structure(Geo,Coord,fs(i),rs(i),Chord(i),I.Material.Wing.UpperPanel,...
              I.Material.Wing.LowerPanel,I.Material.Wing.FrontSpar,I.Material.Wing.RearSpar,I.Structure.Wing.UpperPanelEfficiency,Geo.Wing.RibPitch,Load,y(i));

   end   
   
   for i = 1 :length(y)-1
       
       yst(i) = (y(i)+y(i+1))/2;

% %        if isnan(Au(i))
% % 
% %             for j1=i+1:length(y)
% %                if ~isnan(Au(j1))
% %                    
% %                    tsu(i) = interp1([y(i) y(j1)],[tsu(i-1) tsu(j1)],(y(i)+y(j1))/2,'spline');
% %                    tsl(i) = interp1([y(i) y(j1)],[tsl(i-1) tsl(j1)],(y(i)+y(j1))/2,'spline');
% %                    tfs(i) = interp1([y(i) y(j1)],[tfs(i-1) tfs(j1)],(y(i)+y(j1))/2,'spline');
% %                    trs(i) = interp1([y(i) y(j1)],[trs(i-1) trs(j1)],(y(i)+y(j1))/2,'spline');
% %                            
% %                    Au(i) = Su(i) * tsu(i);
% %                    Al(i) = Sl(i) * tsl(i);
% %                    Afs(i) = h_fs(i) * tfs(i);
% %                    Ars(i) = h_rs(i) * trs(i);
% %                            
% %                    break;
% %                 end
% %            end
% %        end


       if isnan(Au(i)) && sum(isnan(Au))/length(Au) <= 0.3 && i~=1
           nonan = find(~isnan(Au));
           ynonan = y(nonan);
           tsunonan = tsu(nonan);
           tslnonan = tsl(nonan);
           tfsnonan = tfs(nonan);
           trsnonan = trs(nonan);
           
           tsu(i) = interp1(ynonan,tsunonan,y(i),'spline');
           tsl(i) = interp1(ynonan,tslnonan,y(i),'spline');
           tfs(i) = interp1(ynonan,tfsnonan,y(i),'spline');
           trs(i) = interp1(ynonan,trsnonan,y(i),'spline');
           
           Au(i) = Su(i) * tsu(i);
           Al(i) = Sl(i) * tsl(i);
           Afs(i) = h_fs(i) * tfs(i);
           Ars(i) = h_rs(i) * trs(i);
           
       end
       
       Wu(i)  = I.Material.Wing.UpperPanel.rho*g*Au(i)*(y(i+1)-y(i))*Geo.Wing.Span/cos(Geo.Wing.Sweep_half);
       Wl(i)  = I.Material.Wing.LowerPanel.rho*g*Al(i)*(y(i+1)-y(i))*Geo.Wing.Span/cos(Geo.Wing.Sweep_half);
       Wfs(i) = I.Material.Wing.FrontSpar.rho*g*Afs(i)*(y(i+1)-y(i))*Geo.Wing.Span/cos(Geo.Wing.Sweep_half);
       Wrs(i) = I.Material.Wing.RearSpar.rho*g*Ars(i)*(y(i+1)-y(i))*Geo.Wing.Span/cos(Geo.Wing.Sweep_half);
       
       W_box(i) = Wu(i) + Wl(i) + Wfs(i) + Wrs(i);
       clear Choord yu yl
   end 

    W_wingbox = sum(W_box)/9.81;
%     W_w = 7.6621*W_wingbox^0.8509;
    W_w = 10.147*W_wingbox^0.8162;

  
    if I.Disp ==1
        
        disp(' ')
        disp('   y/(b/2)   Chord[m]  tu[mm]     tl[mm]     tfs[mm]   trs[mm]  ')
        disp([yst' Chord tsu.*1000 tsl.*1000 tfs.*1000 trs.*1000])
        disp(' ')
        disp(['Wing Weight   ' num2str(W_w) ' kg']) 
        disp(' ')
    end

    if isnan(W_w)
        W_w = NaN;
        break;
    end
    
    if abs(W_w - W_w0)/W_w < 0.005
        if I.Disp ==1
           disp('Solution is converged');
        end
        break;
    else
        if iter ==20
            W_w = NaN;
        else
            W_w0 = W_w;
        end
    end

end


W_wing.W_Wing = W_w;
W_wing.W_Wingbox = W_wingbox;

W_wingbox
W_w

St.y = yst;
St.tu = tsu;
St.tl = tsl;
St.tfs = tfs;
St.trs = trs;
St.Chord = Chord;

end

%%
function Load = Loads(I,Geo,AS,Wwing,Wfuel)

%% Initializing
% Wwing.y = [0 1];
% Wwing.w = [2*1545/Geo.Wing.Span 0];
% 
% Wfuel = (I.Weight.MTOW - I.Weight.ZFW)*0.5;

Poption = 0;
%%

% q = 0.5*AS.FlightCond.rho*AS.FlightCond.V^2;

%Y = AS.AeroLoads.Wing1.Yst;
%y = linspace(0,1,round((Geo.Wing.Span/2)/Geo.Wing.RibPitch));
y = linspace(0,1,50);


Af      = zeros(length(y),1);
L_aero  = zeros(length(y),1);
L_wing  = zeros(length(y),1);
L_pp    = zeros(length(y),1);
L_fuel  = zeros(length(y),1);


%% fuel tank volume
yf = linspace(I.WingFuelTank.Ystart,I.WingFuelTank.Yend,10);
Chordf = interp1(Geo.Wing.Yle,Geo.Wing.Chord,yf*Geo.Wing.Span/2);
fsf = interp1(Geo.Wing.Yle,Geo.Wing.FS,yf*Geo.Wing.Span/2);
rsf = interp1(Geo.Wing.Yle,Geo.Wing.RS,yf*Geo.Wing.Span/2);
A = zeros(length(yf),1);

for i=1:length(yf)
    
    x = Geo.Wing.Airfoil_coord(1,:,1).*Chordf(i);
    [yu yl] = Airfoil_interp(Geo.Wing.Airfoil_coord,I.Wing(1).AirfoilPosition,yf(i));
    yu = yu.*Chordf(i);
    yl = yl.*Chordf(i);
        
    SPu = spline(x,yu);
    SPl = spline(x,yl);
        
    A(i) = quad(@(x)ppval(SPu,x),fsf(i)*Chordf(i),rsf(i)*Chordf(i)) - quad(@(x)ppval(SPl,x),fsf(i)*Chordf(i),rsf(i)*Chordf(i));
           
end

for i=1:length(y)
    if y(i)>= I.WingFuelTank.Ystart && y(i)< I.WingFuelTank.Yend
        Af(i) = interp1(yf,A,y(i));
    else
        Af(i) = 0;
    end
end


SPa = spline(yf*Geo.Wing.Span/2,A);
Vfuel = quad(@(yf)ppval(SPa,yf),yf(1)*Geo.Wing.Span/2,yf(end)*Geo.Wing.Span/2);
        
% Loads

for i = 1:length(y)
    
    % aerodynamic
    L_aero(i) = interp1(AS.Y, AS.L ,y(i));
    
    % weight
    L_wing(i) = interp1(Wwing.y,-I.n_max*9.81*Wwing.w,y(i));
       
    % power plant
    if I.PP(1).WingEngineNumber~=0 && isfield(I.PP(1),'EnginePosition')
         for j = 1:I.PP(1).WingEngineNumber/2
             if y(i)<=I.PP(j).EnginePosition/(Geo.Wing.Span/2) && y(i+1)>=I.PP(j).EnginePosition/(Geo.Wing.Span/2)
                 L_pp(i)  = - I.n_max*I.PP(j).EngineWeight*9.81;
             end
         end
    end

    % fuel
    
    if y(i)>= I.WingFuelTank.Ystart && y(i)< I.WingFuelTank.Yend
           
        if y(i+1)>I.WingFuelTank.Yend
            y2 = I.WingFuelTank.Yend;
            Af2 = A(end);
        else
            y2 = y(i+1);
            Af2 = Af(i+1);
        end
        
        L_fuel(i) = L_fuel(i+1) + ((Af2+Af(i))/2*(y2-y(i))*Geo.Wing.Span/2)/Vfuel*-I.n_max *9.81*Wfuel;
        
    end
   
end


L = L_aero + L_wing + L_pp + L_fuel;    


T_aero = AS.M;


S_aero = zeros(length(y),1);
S_wing =  zeros(length(y),1);
S_fuel =  zeros(length(y),1);
S_pp =  zeros(length(y),1);


for i = length(y)-1:-1:1
    S_aero(i) = S_aero(i+1) + (L_aero(i)+L_aero(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
    S_wing(i) = S_wing(i+1) + (L_wing(i)+L_wing(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
    S_fuel(i) = S_fuel(i+1) + (L_fuel(i)+L_fuel(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
    S_pp(i)   = S_pp(i+1)   + (L_pp(i)+L_pp(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
end

S = S_aero + S_wing + S_fuel + S_pp;


M_aero = zeros(length(y),1);
M_wing =  zeros(length(y),1);
M_fuel =  zeros(length(y),1);
M_pp =  zeros(length(y),1);

for i = length(y)-1:-1:1
    M_aero(i) = M_aero(i+1) + (S_aero(i)+S_aero(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
    M_wing(i) = M_wing(i+1) + (S_wing(i)+S_wing(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
    M_fuel(i) = M_fuel(i+1) + (S_fuel(i)+S_fuel(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
    M_pp(i)   = M_pp(i+1)   + (S_pp(i)+S_pp(i+1))/2*(y(i+1)-y(i))*Geo.Wing.Span/2;
end

M = M_aero + M_wing + M_fuel + M_pp;
 
T = interp1(AS.Y, T_aero,y);

if Poption == 1
    figure
    hold on
    plot(y,S_aero,'-b');
    plot(y,S_wing,'-r');
    plot(y,S_fuel,'-g');
    plot(y,S_pp,'-c');
    plot(y,S,'-k');
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

end

%%

function [Au Al Afs Ars tsu tsl tfs trs Su Sl h_fs h_rs] = Section_Structure(Geo,Coord,fs,rs,Chord,Mat_u,Mat_l,Mat_fs,Mat_rs,F,RP,Load,ys,Init)

plotoption = 0;

%% Constants

iter_max = 1000;

delta_fs = 1e-6;
delta_rs = 1e-6;
delta_su = 1e-6;
delta_sl = 1e-6;

tsg = 0.0008;

%% airfoil coordinate

x = Coord(:,1);
yu = Coord(:,2);
yl = Coord(:,3);

x = x.*Chord;
yu = yu.*Chord;
yl = yl.*Chord;

Chord0 = Chord/cos(Geo.Wing.Sweep_half);

%% geometry

tmax = max(yu-yl);                                   % maximum thickness

wb = 0;

for i=1:length(x)-1
    
   if x(i)<=rs*Chord && x(i+1)>=fs*Chord
        
        wb = wb+1;
        
        dsu(wb) = sqrt((x(i+1)-x(i))^2+(yu(i+1)-yu(i))^2);
        dx(wb) = (x(i)-x(i+1));
        dsl(wb) = sqrt((x(i+1)-x(i))^2+(yl(i+1)-yl(i))^2);
        
        xx(wb)  = 0.5*(x(i+1)+x(i));
        yxu(wb) = 0.5*(yu(i+1)+yu(i));
        yxl(wb) = 0.5*(yl(i+1)+yl(i));
              
        dydx_u(wb) = (yu(i+1)-yu(i))/(x(i+1)-x(i));
        dydx_l(wb) = (yl(i+1)-yl(i))/(x(i+1)-x(i));
        
        

              
    end
end

    
Su = sum(dsu);
Sl = sum(dsl);

hfs = interp1(x,yu,fs*Chord) - interp1(x,yl,fs*Chord);
hrs = interp1(x,yu,rs*Chord) - interp1(x,yl,rs*Chord);


yfs = linspace(interp1(x,yl,fs*Chord),interp1(x,yu,fs*Chord),21);
yrs = linspace(interp1(x,yl,rs*Chord),interp1(x,yu,rs*Chord),21);

yxfs = zeros(20,1);
yxrs = zeros(20,1);
dsfs = zeros(20,1);
dsrs = zeros(20,1);

for i= 1:20
    yxfs(21-i) = (yfs(i+1)+yfs(i))/2;
    yxrs(i) = (yrs(i+1)+yrs(i))/2;
    
    dsfs(21-i) = yfs(i+1)-yfs(i);
    dsrs(i) = yrs(i+1)-yrs(i);
end
    


xfs = fs*Chord*ones(20,1);
xrs = rs*Chord*ones(20,1);

%% Initializing

if nargin < 14
   
    eta = 0.8;
    tsu = 0.025*tmax;
    tsl = 0.025*tmax;
    tfs = 0.025*tmax;
    trs = 0.025*tmax;
else
    eta = Init.eta;
    tsu = Init.tsu;
    tsl = Init.tsl;
    tfs = Init.tfs;
    trs = Init.trs;
end
    


%% Iteration for skin thickness
Conv = 0;

SPu = spline(x,yu);
SPl = spline(x,yl);
        
A_box = quad(@(x)ppval(SPu,x),fs*Chord,rs*Chord) - quad(@(x)ppval(SPl,x),fs*Chord,rs*Chord);

for iter = 1:iter_max

    yun  = yxu - tsu/2;
    yln  = yxl + tsl/2;
    
    xfs = xfs + tfs/2;
    xrs = xrs - trs/2;


   %% Shear and Bending moment centre
   
    xxl = zeros(1,length(xx));
    for i=1:length(xx)
        xxl(i) = xx(end+1-i);
    end
    
    dsq = [dsu dsfs' dsl dsrs'];
    yq = [yun yxfs' yln yxrs'];
    xq = [xx xfs' xxl xrs'];
    tsq = [tsu*ones(1,length(dsu)) tfs*ones(1,length(dsfs)) tsl*ones(1,length(dsl)) trs*ones(1,length(dsrs))];
    
    Cu = yun - dydx_u.*xx;
    Cl = yln - dydx_l.*xxl;
    pu = abs(Cu./sqrt(1+dydx_u.^2));
    pl = abs(Cl./sqrt(1+dydx_l.^2));
    pfs = xfs;
    prs = xrs;
    
    P = [pu -pfs' pl prs'];

    Ixx = sum((dsq.*tsq).*yq.^2);
    Iyy = sum((dsq.*tsq).*xq.^2);
    Ixy = sum((dsq.*tsq).*(xq.*yq));
    
    qbs = zeros(1,length(dsq));
    qbs(1) = Ixy/(Ixx*Iyy-Ixy^2)*tsq(1)*xq(1)*dsq(1) - Iyy/(Ixx*Iyy-Ixy^2)* tsq(1)*yq(1)*dsq(1);
    
    for s=2:length(dsq)
         qbs(s) = qbs(s-1) + Ixy/(Ixx*Iyy-Ixy^2)*tsq(s)*xq(s)*dsq(s) - Iyy/(Ixx*Iyy-Ixy^2)* tsq(s)*yq(s)*dsq(s);
    end
     
   
    qs0 = -sum((qbs./tsq).*dsq)/sum(dsq./tsq);
    
    qs = qbs+qs0;
    xs = sum((qs.*P).*dsq);
    
    if xs/Chord >0.6 || xs/Chord <0.2
        xs = 0.40*Chord;
    end
    
    y0   = (sum((yun*tsu).*dsu)+sum((yln*tsl).*dsl)+sum((yxfs*tfs).*dsfs)+sum((yxrs*trs).*dsrs))...
        /(sum(tsu*dsu)+sum(tsl*dsl)+sum(tfs*dsfs)+sum(trs*dsrs));

    
    
    %% Loads
    
    
    Sweep = atan((Geo.Wing.Xle(end)+xs/Chord*Geo.Wing.Chord(end) - Geo.Wing.Xle(1)+xs/Chord*Geo.Wing.Chord(1))...
        /(Geo.Wing.Span/2));
    
    Chords = Chord0 *cos(Sweep);
    
    S = interp1(Load.y,Load.S,ys)*1.5;
    M = interp1(Load.y,Load.M,ys)/cos(Sweep)*1.5;

    yt = linspace(ys,1,10);
    Tm1 = (interp1(Load.y,Load.L.L,yt).*interp1(Geo.Wing.Yle/(Geo.Wing.Span/2),Geo.Wing.Chord,yt)).*(xs/Chord-0.25);
    
    
%     Tm2 = -(interp1(Load.y,Load.L.L,yt).*(yt-ys)).*tan(Sweep);
%     Tm3 = interp1(Load.y,Load.T,yt);
%     Tm = Tm1 + Tm2 + Tm3;
%     
%     T = 0;
%     for m=9:-1:1
%         T = T + (Tm(m+1)+Tm(m))/2*(yt(m+1)-yt(m))*Geo.Wing.Span/2;
%     end
    %T =T/cos(Sweep)*1.5;
        
 T = (interp1(Load.y,Load.T,ys)+Tm1)/cos(Sweep)*1.5;
  

        
    %% Panel thicknesses
    
    ymaxu = max(abs(yun - y0));
    ymaxl = max(abs(yln-y0));
           
    sigma_max_b = F*(M*Mat_u.E/(eta*tmax*(rs-fs)*Chords*RP))^0.5;
    sigma_maxl = Mat_l.Sigma_tensile;
    sigma_max_c = Mat_u.Sigma_compressive;
    sigma_maxu  = min(sigma_max_c,sigma_max_b);
    
    
    eta = 1/(tmax)*(sum(((yun-y0).^2).*dsu/(Su*ymaxu)) + sum(((yln-y0).^2).*dsl/(Sl*ymaxl)));
     
     
    tsu_n = max([(M/sigma_maxu)/(eta*tmax*Su) tsg]);
    tsl_n = max([(M/sigma_maxl)/(eta*tmax*Sl) tsg]);
    
    ru = abs(tsu - tsu_n); 
    rl = abs(tsl - tsl_n);

    tsu = tsu_n;
    tsl = tsl_n;
     
    dAu = dsu.*tsu;
    dAl = dsl.*tsl;

    Au = sum(dAu);
    Al = sum(dAl);
    
    %% Spar thicknesses

    Taw_fs = Mat_fs.Sigma_compressive/2;
    Taw_rs = Mat_rs.Sigma_compressive/2;
    
    % shear flows
    
    h_fs = hfs - tsu - tsl;
    h_rs = hrs - tsu - tsl;
    
    Ssfs = h_fs^2/(h_fs^2+h_rs^2)*S;
    Ssrs = h_rs^2/(h_fs^2+h_rs^2)*S;
    
    Qs_fs = Ssfs/h_fs;
    Qs_rs = Ssrs/h_rs;
    
    Qt = T/(2*A_box);
    
    Qfs = Qs_fs + Qt;
    Qrs = Qs_rs - Qt;
    
    tfs_n = max([abs(Qfs/Taw_fs) tsg]);
    trs_n = max([abs(Qrs/Taw_rs) tsg]);

    Afs = h_fs * tfs_n;
    Ars = h_rs * trs_n;

    rfs = abs(tfs_n - tfs);
    rrs = abs(trs_n - trs);
    
    tfs = tfs_n;
    trs = trs_n;

    %% convergency?
    
    if ru < delta_su && rl < delta_sl && rfs < delta_fs && rrs < delta_rs && ~isnan(xs)
        Conv = 1;
        break;
    else
%        clear y0 x0 dAu dAu Au Al yun yln ymaxu ymaxl ymax sigma_maxu
%        clear Afs Ars tfs_n trs_n Sfs Srs
    end
    
end


if Conv ==0
    Au  = NaN;
    Al  = NaN;
    tsu = NaN;
    tsl = NaN;
    
    Afs = NaN;
    Ars = NaN;
    tfs = NaN;
    trs = NaN;
end


end
 


  
            
            
        




    




    


    


