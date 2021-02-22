function [Cleff, Cdpeff, Cdfeff, dCleff_daeff, dCleff_dMeff, dCleff_dReeff, dCdpeff_daeff, dCdpeff_dMeff, ...
    dCdpeff_dReeff, dCdfeff_daeff, dCdfeff_dMeff, dCdfeff_dReeff, dCleff_dX, dCdpeff_dX, dCdfeff_dX] = ...
    Airfoil_analysis_parallel(AC,n2d,Eff,CSTi,Mode,iter)

global Restart

ModeNum  = AC.Wing.Airfoils.ModeNum;

aeff = Eff(1:n2d);
Meff = Eff(n2d+1:2*n2d);
Reeff = Eff(2*n2d+1:3*n2d);

Mode = reshape(Mode',2*ModeNum,n2d)';

%% MSES analysis

Tran = 'ForceT';
ModeNum = size(AC.Wing.Airfoils.Chebyshev,2)/2;

parfor i=1:n2d  
    Flow{i}.M    = Meff(i);
    Flow{i}.Re   = Reeff(i); 
    Flow{i}.Alfa = aeff(i)*180/pi;
    Flow{i}.CL   = 0; 
    Flow{i}.Mode = 5;  
end


if iter <= Restart
     ! rm -r Tstorage*
    parfor i=1:n2d 
        [~,~] = system(['cp -avr Storage Tstorage' num2str(i)]); 
    end

    parfor i=1:n2d 
       Res{i} = mses_analysis(Flow{i},ModeNum,Tran,i,Mode(i,:),iter,CSTi(i,:)); 
    end

else
    parfor i=1:n2d  
        
        cd(['Tstorage' num2str(i)])
        fid=fopen('params.test','w+');

        fprintf(fid,'%g 0\n',2*ModeNum);
        for j=1:2*ModeNum
            fprintf(fid,'%g \n', Mode(i,j));
        end

        fprintf(fid,'%g %g %g\n', round(aeff(i)*180/pi*1000000)/1000000, 0,round(Meff(i)*1000000)/1000000); %round(cl2d(i)*100)/100
        fprintf(fid,'%g \n', round(Reeff(i)/1e7*1000000)*1e7/1000000);

%         fprintf(fid,'%6.10f %6.10f %6.10f\n', aeff(i)*180/pi, 0,Meff(i)); %round(cl2d(i)*100)/100
%         fprintf(fid,'%6.10f \n', Reeff(i));

        fclose(fid);
        cd ..

        R = MSES_Run(1,ModeNum,i);
        
        if ~isnan(R.M)
            Res{i} = R;
        else  
            [~,~] = system(['rm -r Tstorage' num2str(i)]); 
            [~,~] = system(['cp -avr Storage Tstorage' num2str(i)]); 
            
            Res{i} = mses_analysis(Flow{i},ModeNum,Tran,i,Mode(i,:),iter,CSTi(i,:)); 
        end
        
 
    end
    
end


%%
% Res.Sens: {2}=CL, {3}=CM, {4}=CD, {5}=CDw, {6}=CDv, {7}=CDf

Cleff = zeros(n2d,1);
Cdpeff = zeros(n2d,1);
Cdfeff = zeros(n2d,1);

dCleff_daeff = zeros(n2d);   dCleff_dMeff = zeros(n2d);   dCleff_dReeff = zeros(n2d);
dCdpeff_daeff = zeros(n2d);  dCdpeff_dMeff = zeros(n2d);  dCdpeff_dReeff = zeros(n2d);
dCdfeff_daeff = zeros(n2d);  dCdfeff_dMeff = zeros(n2d);  dCdfeff_dReeff = zeros(n2d);

dCleff_dX = zeros(n2d,n2d*2*ModeNum);  dCdpeff_dX = zeros(n2d,n2d*2*ModeNum);   dCdfeff_dX = zeros(n2d,n2d*2*ModeNum);

for i=1:n2d
    R = Res{i};

    Cleff(i) = R.CL;
    Cdpeff(i) = R.CDp;
    Cdfeff(i) = R.CDf;

    if ~isempty(R.Sens.Alfa{2})
        dCleff_daeff(i,i)  = R.Sens.Alfa{2}*180/pi;
        dCdpeff_daeff(i,i) = (R.Sens.Alfa{4}-R.Sens.Alfa{7})*180/pi;
        dCdfeff_daeff(i,i) = R.Sens.Alfa{7}*180/pi;

        dCleff_dMeff(i,i)  =  R.Sens.Mach{2};
        dCdpeff_dMeff(i,i) =  R.Sens.Mach{4}-R.Sens.Mach{7};
        dCdfeff_dMeff(i,i) =  R.Sens.Mach{7};

        dCleff_dReeff(i,i) = R.Sens.Re{2}/Reeff(i);
        dCdpeff_dReeff(i,i) = R.Sens.Re{4}/Reeff(i) - R.Sens.Re{7}/Reeff(i);
        dCdfeff_dReeff(i,i) = R.Sens.Re{7}/Reeff(i);

        for m = 1:2*ModeNum
            dCleff_dX(i, (i-1)*2*ModeNum + m)  = R.Sens.Mode{m}{2};
    %         dCdpeff_dX(i, (i-1)*2*ModeNum + m) = R.Sens.Mode{m}{4}-R.Sens.Mode{m}{7};
    %         dCdfeff_dX(i, (i-1)*2*ModeNum + m) = R.Sens.Mode{m}{7};

            % **** Approximating sensitivities because MSES does not report dCDf/dMode correctly
            dCdpeff_dX(i, (i-1)*2*ModeNum + m) = (Cdpeff(i)/(Cdfeff(i) + Cdpeff(i)))*R.Sens.Mode{m}{4};
            dCdfeff_dX(i, (i-1)*2*ModeNum + m) = (Cdfeff(i)/(Cdfeff(i) + Cdpeff(i)))*R.Sens.Mode{m}{4};

        end
    end

    
end
end