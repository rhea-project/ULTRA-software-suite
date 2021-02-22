function Res = mses_analysis(Flow,ModeNum,Tran,SecNum,Mode,iter,varargin)

Res=[];

% if nargin ==5 || nargin==6;
%     element =1;
% elseif nargin ==8
%     element = 2;
% elseif nargin ==10
%     element =3;
% end

element = 1;

cd(['Tstorage' num2str(SecNum)])

%% Grid parameters

if element ==1
    [~,~] = system('cp gridpar_single.test gridpar.test');
elseif element ==2
    [~,~] = system('cp gridpar_flap.test gridpar.test');
end


%% Create blade.air file

fid = fopen('blade.test','w+');

fprintf(fid,'Airfoil\n');
fprintf(fid,'-2.0 3.0 -3.0 3.5\n');


if nargin ==7
    Nx = 100;
    x = [0 1-cos(pi/Nx:pi/Nx:pi)]/2;

    A = varargin{1};
    
    Au = A(1:length(A)/2);
    Al = A(length(A)/2+1:end);

    yu = CSTairfoil(Au,x,0.00);
    yl = CSTairfoil(Al,x,-0.00);

    X = [flipud(x'); x(2:end)'];
    Y = [flipud(yu); yl(2:end)];
    
    for i=1:length(X)
        fprintf(fid,'%g %g\n',X(i),Y(i));
    end
    fprintf(fid,'999 999\n');
    
else

    for e=1:element
        X = varargin{1+(e-1)*2};
        Y = varargin{2+(e-1)*2};
        for i=1:length(X)
            fprintf(fid,'%g %g\n',X(i),Y(i));
        end
        fprintf(fid,'999 999\n');
    end

end

fclose(fid);

    

%% Create  mses file

fid=fopen('mses.test','w+');

fprintf(fid,'3 4 5 7 10 15 20\n');
if Flow.Mode ~= 0
    fprintf(fid,'3 4 %g 7 15 17 20\n', Flow.Mode); % 5 for AOA 6 for CL
else
    fprintf(fid,'3 4 5 7 15 17 20\n');
end    
fprintf(fid,'%g %g %g\n',round(Flow.M*1000000)/1000000,round(Flow.CL*1000000)/1000000,round(Flow.Alfa*1000000)/1000000);
fprintf(fid,'4 2\n');
fprintf(fid,'%g 9.0\n',round(Flow.Re/1e7*1000000)*1e7/1000000);
if strcmpi(Tran,'ForceT')
    fprintf(fid,'0.02 0.02 \n');  % only for single element airfoil
elseif strcmpi(Tran,'FreeT')
    fprintf(fid,'1 1\n');
end    
fprintf(fid,'0.97  1.2 \n');

fprintf(fid,'1 1\n');
fprintf(fid,'%f 0\n',2*ModeNum);

fclose(fid);

%% Create mode file

fid=fopen('modes.test','w+');

for i= 1:ModeNum
    fprintf(fid,'%g %g 1 0 1 1 \n',i,20+i);
end
for j=1:ModeNum
    fprintf(fid,'%g %g 1 0 -1 1 \n',j+i,20+j);
end

fclose(fid);

%% Create Param file

fid=fopen('params.test','w+');

fprintf(fid,'%g 0\n',2*ModeNum);
for j=1:2*ModeNum
    fprintf(fid,'%g \n', Mode(j));
end

fprintf(fid,'%g %g %g\n', round(Flow.Alfa*1000000)/1000000, Flow.CL,round(Flow.M*1000000)/1000000); %round(cl2d(i)*100)/100
fprintf(fid,'%g \n', round(Flow.Re/1e7*1000000)*1e7/1000000);

% fprintf(fid,'%6.10f %6.10f %6.10f\n', Flow.Alfa, Flow.CL,Flow.M); %round(cl2d(i)*100)/100
% fprintf(fid,'%6.10f \n', Flow.Re);

fclose(fid);

%% Create Mset logfile

fid = fopen('msetlog.dat','w+');

fprintf(fid,'1\n');
fprintf(fid,'999\n'); %,round(Flow.Alfa*100)/100);
fprintf(fid,'2\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
if element ==3
    fprintf(fid,'\n');
end
fprintf(fid,'3\n');
fprintf(fid,'4\n');
fprintf(fid,'0\n');

fclose(fid);

%% Create mplot log file

fid = fopen('mplotlog.dat','w+');

fprintf(fid,'1\n');
fprintf(fid,'12\n');
fprintf(fid,'\n');
% fprintf(fid,'1\n');
% fprintf(fid,'2\n');
% fprintf(fid,'\n');
% fprintf(fid,'1\n');
% fprintf(fid,'6\n');
fprintf(fid,'\n');
fprintf(fid,'0\n');

fclose(fid);  

%% lindop logfile

fid = fopen('lindoplog.dat','w+');

fprintf(fid,'5\n');
fprintf(fid,'a\n');
fprintf(fid,'m\n');
fprintf(fid,'r\n');
for i=1:2*ModeNum
    fprintf(fid,'g%f\n',i);
end
fprintf(fid,'\n');
fprintf(fid,'0\n');

fclose(fid);  


%% creat mpolar

if Flow.Mode ==0
    s = [0 1 2 3 4 4.5 5 5.2 5.4 5.6 5.8 6];       
%     s = [4 4.5 5 5.5 6];
    fid = fopen('spec.test','w+');
    fprintf(fid,'5\n');
    for i=1:length(s)
        fprintf(fid,'%g\n', s(i));
    end
    for i=1:100
        fprintf(fid,'%g\n', s(end)+i/10);
    end
    fclose(fid);
end


%% Run Mses


[~,~] = system('./mset_noplot test < msetlog.dat');
if Flow.Mode ~=0
    [~,~] = system('./mses test -100');  %[~,~] = system('mses test -100');
    [~,~] = system('./mses test +1000 > msesconvergence.dat');  %[~,~] = system('mses test +1000 > msesconvergence.dat');
    
    fid  = fopen('msesconvergence.dat','r');
    temp = textscan(fid,' %s');
    fclose(fid);
    
    if ~strcmpi( temp{1,1}{end-107-(2*ModeNum-20)} , 'Converged' )
        [~,~] = system('rm msesconvergence.dat');
        [~,~] = system('./mset_noplot test < msetlog.dat');
        [~,~] = system('./mses test +1000 > msesconvergence.dat'); %[~,~] = system('mses test +1000 > msesconvergence.dat');
        
        fid  = fopen('msesconvergence.dat','r');
        temp = textscan(fid,' %s');
        fclose(fid);
    end
    
    if ~strcmpi( temp{1,1}{end-107-(2*ModeNum-20)} , 'Converged' )
        Res.M = NaN;
        Res.AOA =  NaN;
        Res.CL = NaN;
        Res.Re = NaN;
        Res.CD =  NaN;
        Res.CM =  NaN;
        Res.LD =  NaN;
        Res.CDv =  NaN;
        Res.CDw =  NaN;
        Res.CDf =  NaN;
        Res.CDp =  NaN;
        Res.CL_alfa =  NaN;
        Res.CD_lafa =  NaN;
        Res.CM_alfa =  NaN;
        Res.dCL_dM = NaN;
        Res.dCD_dM  =  NaN;
        Res.dCM_dM =  NaN;
    else

    system('./mplot test < mplotlog.dat > msesdata.dat'); %system('mplot test < mplotlog.dat > msesdata.dat');
    
        fid = fopen('msesdata.dat','r');
        Res_MSES =textscan(fid,'%s','headerlines',0);
        fclose(fid);

        R = Res_MSES{1};

        for i = 1:length(R)
            if strcmp(R(i),'Ma')
                break;
            end
        end

        i = i-1;

        Res.M = str2double(R(i+3));
        Res.Re = str2double(R(i+10));
        Res.AOA =  str2double(R(i+6));
        Res.CL = str2double(R(i+16));
        Res.CD =  str2double(R(i+19));
        Res.CM =  str2double(R(i+22));
        Res.LD =  str2double(R(i+25));
        Res.CDv =  str2double(R(i+33));
        Res.CDw =  str2double(R(i+36));
        Res.CDf =  str2double(R(i+40));
        Res.CDp =  str2double(R(i+43));
        Res.CL_alfa =  str2double(R(i+63));
        Res.CD_lafa =  str2double(R(i+66));
        Res.CM_alfa =  str2double(R(i+69));
        Res.dCL_dM =  str2double(R(i+72));
        Res.dCD_dM  =  str2double(R(i+75));
        Res.dCM_dM =  str2double(R(i+78));

        system('rm msesdata.dat');

    end
    
else
    system('./mpolar test'); %system('mpolar test');
    
    fid = fopen('polar.test','r');
    Res_Mpolar =textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',13);
    fclose(fid);
    
    alpha = Res_Mpolar{1};
    Cl    = Res_Mpolar{2};
    Cd    = Res_Mpolar{3};
    
    x = linspace(min(alpha),max(alpha),100);
    P = spline(alpha,Cl);
    [CLmax, na] = max(ppval(x,P));
    alpha_max = x(na);
    
    Res.Alpha = alpha;
    Res.CL = Cl;
    Res.CD = Cd;
    Res.CLmax = CLmax;
    Res.Aclmax = alpha_max;
    
    system('rm spec.test');
    system('rm polar.test');
    system('rm polarx.test');
end

%% Reaad sensitivity
system('./lindop test <lindoplog.dat >lindopdata.dat');

fid = fopen('lindopdata.dat','r');
Res.Sens.Alfa =textscan(fid,'%f %f %f %f %f %f %f','headerlines',68);

if ~isempty(Res.Sens.Alfa{1})
    Res.Sens.Mach =textscan(fid,'%f %f %f %f %f %f %f','headerlines',11);
    Res.Sens.Re =textscan(fid,'%f %f %f %f %f %f %f','headerlines',11);
    for i=1:2*ModeNum
        Res.Sens.Mode{i} =textscan(fid,'%f %f %f %f %f %f %f','headerlines',11);
    end
else
    fclose(fid);
    fid = fopen('lindopdata.dat','r');
    Res.Sens.Alfa =textscan(fid,'%f %f %f %f %f %f %f','headerlines',72);
    Res.Sens.Mach =textscan(fid,'%f %f %f %f %f %f %f','headerlines',11);
    Res.Sens.Re =textscan(fid,'%f %f %f %f %f %f %f','headerlines',11);
    for i=1:2*ModeNum
        Res.Sens.Mode{i} =textscan(fid,'%f %f %f %f %f %f %f','headerlines',11);
    end
end
fclose(fid);
    

        
%% removing files

% system('rm mplotlog.dat');
% system('rm msetlog.dat');
system('rm msesconvergence.dat');
% system('rm lindoplog.dat');
system('rm lindopdata.dat');

% system('rm blade.test');
% system('rm mses.test');
% system('rm mdat.test');
% system('rm gridpar.test');



 cd .. 
 
end

function [Y, Rle, beta]=CSTairfoil(A,x,zeta)

N1 = 0.5;   %Class function N1
N2 = 1;     %Class function N2

if nargin<3
    zeta = 0.000;     % TE gap
end
    
n = length(A)-1;

%evaluate required functions for each X-ordinate
for i = 1:length(x)
    
    %calculate Class Function for x(i):
    C(i) = (x(i)^N1)*(1-x(i))^N2;
    
    %calculate Shape Functions surface at x(i)
    S(i) = 0;  %Shape function initially zero
    for j = 0:n
        Krnu = factorial(n)/(factorial(j)*factorial(n-j));
        S(i) = S(i) + A(j+1)*Krnu*(1-x(i))^(n-j)*x(i)^(j);
    end

    if i==1
        Rle=(S(i)^2)/2;
    end
    if i==length(x)
        beta=atan(S(i));
    end
    %calculate upper and lower surface ordinates at x(i)
    Y(i) = C(i)*S(i) + x(i)*zeta;

 

end

Y= Y';

end
