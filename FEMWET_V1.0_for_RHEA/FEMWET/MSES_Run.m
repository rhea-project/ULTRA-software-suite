function Res = MSES_Run(Case,ModeNum,SecNum)

cd(['Tstorage' num2str(SecNum)]);


[~,~] = system('./mses test +1000 > msesconvergence.dat'); %system('mses test +1000 > msesconvergence.dat')

fid  = fopen('msesconvergence.dat','r');
temp = textscan(fid,' %s');
fclose(fid);

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

    system('./mplot test < mplotlog.dat > msesdata.dat'); % system('mplot test < mplotlog.dat > msesdata.dat');

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
system('rm msesconvergence.dat');

if Case ==1
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
    system('rm lindopdata.dat');
end

cd ..
