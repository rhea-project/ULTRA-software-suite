function [Failure, Loc, Z, dg_dx] = Stress(AC,U,X,cNode,Xea,Yea,dcNode,dXea,dYea,DV,DP)

global Adjoint

plotoption =0;


ModeNum  = AC.Wing.Airfoils.ModeNum;
nT       = AC.Structure.nT;
yT       = AC.Structure.yT;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
yAirfoil = AC.Wing.Airfoils.yAirfoil;
nNodea   = AC.FEM.nNodea;
nElm     = AC.FEM.nElm;
Conc     = AC.FEM.Conc;
nSkin    = AC.FEM.nSkin;
nSpar    = AC.FEM.nSpar;
ActiveDOFs = AC.FEM.ActiveDOFs;
nChord   = AC.Q3D.Grid.nChord;




for i=1:nElm
        
    % Displacement vector in local ccordinates
    
    Ue = U((i-1)*6+1:(i-1)*6+12); % element dispacement vector
        
     %  Coord transformation matrix struct(x,y,z) to real(y,z,x)
    struct2real = zeros(6,6);
    struct2real(1,3) = 1;
    struct2real(2,1) = 1;
    struct2real(3,2) = 1;
    struct2real(4,6) = 1;
    struct2real(5,4) = 1;
    struct2real(6,5) = 1;
    s2r = zeros(12,12);
    s2r(1:6, 1:6) = struct2real;
    s2r(7:12, 7:12) = struct2real;

    Ue = s2r \ Ue;
    
    %%
    
    if Adjoint ==1
        [g_tension_u(i,:), g_compression_u(i,:), g_buckling_u(i,:), g_tension_l(i,:), g_compression_l(i,:), g_buckling_l(i,:), g_shear_fs(i,:), g_shear_rs(i,:), g_buckling_fs(i,:),g_buckling_rs(i,:), ... 
            Au(i), Al(i), Afs(i), Ars(i), Xu(i,:), Yu(i,:), Xl(i,:), Yl(i,:), Xfs(i,:), Yfs(i,:), Xrs(i,:), Yrs(i,:),  ...
             dg_tension_du_u, dg_tension_du_l, dg_compression_du_u, dg_compression_du_l, dg_buckling_du_u, dg_buckling_du_l, dg_shear_du_fs, dg_shear_du_rs, dg_buckling_du_fs, dg_buckling_du_rs, ...
             dg_tension_dx_u, dg_tension_dx_l, dg_compression_dx_u, dg_compression_dx_l, dg_buckling_dx_u, dg_buckling_dx_l, dg_shear_dx_fs, dg_shear_dx_rs, dg_buckling_dx_fs, dg_buckling_dx_rs, dAu_dx, dAl_dx, dAfs_dx, dArs_dx] ...
             = Element_Stress(AC,Ue,nNodea,Conc,i,X,nT,yT,nAirfoil,yAirfoil,nSkin,nSpar,s2r,cNode,Xea,Yea,dcNode,dXea,dYea,DV,DP);
    else
        [g_tension_u(i,:), g_compression_u(i,:), g_buckling_u(i,:), g_tension_l(i,:), g_compression_l(i,:), g_buckling_l(i,:), g_shear_fs(i,:), g_shear_rs(i,:), g_buckling_fs(i,:),g_buckling_rs(i,:), ...
             Au(i), Al(i), Afs(i), Ars(i), Xu(i,:), Yu(i,:), Xl(i,:), Yl(i,:), Xfs(i,:), Yfs(i,:), Xrs(i,:), Yrs(i,:)] = Element_Stress(AC,Ue,nNodea,Conc,i,X,nT,yT,nAirfoil,yAirfoil,nSkin,nSpar,s2r,cNode,Xea,Yea,dcNode,dXea,dYea,DV,DP);
    end

%      Xu(i,:) = Xu(i,:)   + (cNode(i,1)+cNode(i+1,1))/2;
%      Yu(i,:) = Yu(i,:)   + (cNode(i,3)+cNode(i+1,3))/2;
%      Xl(i,:) = Xl(i,:)   + (cNode(i,1)+cNode(i+1,1))/2;
%      Yl(i,:) = Yl(i,:)   + (cNode(i,3)+cNode(i+1,3))/2;
%      Xfs(i,:) = Xfs(i,:) + (cNode(i,1)+cNode(i+1,1))/2;
%      Yfs(i,:) = Yfs(i,:) + (cNode(i,3)+cNode(i+1,3))/2;
%      Xrs(i,:) = Xrs(i,:) + (cNode(i,1)+cNode(i+1,1))/2;
%      Yrs(i,:) = Yrs(i,:) + (cNode(i,3)+cNode(i+1,3))/2;
     

%% Adjoint

    if Adjoint ==1

        % ***** lambda vector *****
        
        z_tension_u     = zeros(length(U),nSkin);
        z_compression_u = zeros(length(U),nSkin);
        z_buckling_u    = zeros(length(U),nSkin);
       
        z_tension_l     = zeros(length(U),nSkin);
        z_compression_l = zeros(length(U),nSkin);
        z_buckling_l    = zeros(length(U),nSkin);       
        
        z_shear_fs      = zeros(length(U),nSpar);
        z_buckling_fs   = zeros(length(U),nSpar);
        
        z_shear_rs      = zeros(length(U),nSpar);
        z_buckling_rs   = zeros(length(U),nSpar);
               
        
        z_tension_u((i-1)*6+1:(i-1)*6+12,:)     = dg_tension_du_u;
        z_compression_u((i-1)*6+1:(i-1)*6+12,:) = dg_compression_du_u;
        z_buckling_u((i-1)*6+1:(i-1)*6+12,:)    = dg_buckling_du_u;
        
        z_tension_l((i-1)*6+1:(i-1)*6+12,:)     = dg_tension_du_l;
        z_compression_l((i-1)*6+1:(i-1)*6+12,:) = dg_compression_du_l;
        z_buckling_l((i-1)*6+1:(i-1)*6+12,:)    = dg_buckling_du_l;
        
        z_shear_fs((i-1)*6+1:(i-1)*6+12,:)      = dg_shear_du_fs;
        z_shear_rs((i-1)*6+1:(i-1)*6+12,:)      = dg_shear_du_rs;
        
        z_buckling_fs((i-1)*6+1:(i-1)*6+12,:)   = dg_buckling_du_fs;
        z_buckling_rs((i-1)*6+1:(i-1)*6+12,:)   = dg_buckling_du_rs;
        
               
        Z.U.tension(i,:,:)     = [zeros(nElm*nChord,nSkin) ; z_tension_u(ActiveDOFs,:) ; zeros(1,nSkin)];
        Z.U.compression(i,:,:) = [zeros(nElm*nChord,nSkin) ; z_compression_u(ActiveDOFs,:) ; zeros(1,nSkin)];
        Z.U.buckling(i,:,:)    = [zeros(nElm*nChord,nSkin) ; z_buckling_u(ActiveDOFs,:) ; zeros(1,nSkin)];

        Z.L.tension(i,:,:)     = [zeros(nElm*nChord,nSkin) ; z_tension_l(ActiveDOFs,:) ; zeros(1,nSkin)];
        Z.L.compression(i,:,:) = [zeros(nElm*nChord,nSkin) ; z_compression_l(ActiveDOFs,:) ; zeros(1,nSkin)];
        Z.L.buckling(i,:,:)    = [zeros(nElm*nChord,nSkin) ; z_buckling_l(ActiveDOFs,:) ; zeros(1,nSkin)];

        Z.FS.shear(i,:,:)      = [zeros(nElm*nChord,nSpar) ; z_shear_fs(ActiveDOFs,:) ; zeros(1,nSpar)];
        Z.FS.buckling(i,:,:)   = [zeros(nElm*nChord,nSpar) ; z_buckling_fs(ActiveDOFs,:) ; zeros(1,nSpar)];

        Z.RS.shear(i,:,:)      = [zeros(nElm*nChord,nSpar) ; z_shear_rs(ActiveDOFs,:) ; zeros(1,nSpar)];
        Z.RS.buckling(i,:,:)   = [zeros(nElm*nChord,nSpar) ; z_buckling_rs(ActiveDOFs,:) ; zeros(1,nSpar)];

       
        % *** dg_dx *****
        
        dg_dx.U.tension(i,:,:)       = dg_tension_dx_u;
        dg_dx.U.compression(i,:,:)   = dg_compression_dx_u;
        dg_dx.U.buckling(i,:,:)      = dg_buckling_dx_u;
        
        dg_dx.L.tension(i,:,:)       = dg_tension_dx_l;
        dg_dx.L.compression(i,:,:)   = dg_compression_dx_l;
        dg_dx.L.buckling(i,:,:)      = dg_buckling_dx_l;
        
        dg_dx.FS.shear(i,:,:)        = dg_shear_dx_fs;
        dg_dx.FS.buckling(i,:,:)     = dg_buckling_dx_fs;
        
        dg_dx.RS.shear(i,:,:)        = dg_shear_dx_rs;
        dg_dx.RS.buckling(i,:,:)     = dg_buckling_dx_rs;
        

        % ***** dW_dx *****
        
%         dWu_dx(i,:)  = dAu_dx' * AC.St.Mat_u.rho*(yNode(i+1)-yNode(i))*b/cos(Sweep_half) + ...
%                         Au(i)*AC.St.Mat_u.rho*(yNode(i+1)-yNode(i))/cos(Sweep_half)*db_dx + ...
%                         Au(i)*AC.St.Mat_u.rho*(yNode(i+1)-yNode(i))*b*sin(Sweep_half)/(cos(Sweep_half)^2)*dSweep_dx; 
%                     
%         dWl_dx(i,:)  = dAl_dx' * AC.St.Mat_l.rho*(yNode(i+1)-yNode(i))*b/cos(Sweep_half) + ...
%                         Al(i)*AC.St.Mat_l.rho*(yNode(i+1)-yNode(i))/cos(Sweep_half)*db_dx + ...
%                         Al(i)*AC.St.Mat_l.rho*(yNode(i+1)-yNode(i))*b*sin(Sweep_half)/(cos(Sweep_half)^2)*dSweep_dx; 
%  
%         dWfs_dx(i,:) = dAfs_dx' * AC.St.Mat_fs.rho*(yNode(i+1)-yNode(i))*b/cos(Sweep_half) + ...
%                         Afs(i)*AC.St.Mat_fs.rho*(yNode(i+1)-yNode(i))/cos(Sweep_half)*db_dx + ...
%                         Afs(i)*AC.St.Mat_fs.rho*(yNode(i+1)-yNode(i))*b*sin(Sweep_half)/(cos(Sweep_half)^2)*dSweep_dx; 
%         
%         dWrs_dx(i,:) = dArs_dx' * AC.St.Mat_rs.rho*(yNode(i+1)-yNode(i))*b/cos(Sweep_half) + ...
%                         Ars(i)*AC.St.Mat_rs.rho*(yNode(i+1)-yNode(i))/cos(Sweep_half)*db_dx + ...
%                         Ars(i)*AC.St.Mat_rs.rho*(yNode(i+1)-yNode(i))*b*sin(Sweep_half)/(cos(Sweep_half)^2)*dSweep_dx; 
% 
%         dWt_dx(i,:)  = dWu_dx(i,:) + dWl_dx(i,:) + dWfs_dx(i,:) + dWrs_dx(i,:);
          dWt_dx = 0;  
 
    else
        lambda = NaN;
        dg_dx = NaN;
    end
end


%% Output

% Loc.Y = cNode(1:end-1,2);
% Loc.Xu = Xu;
% Loc.Zu = Yu;
% Loc.Xl = Xl;
% Loc.Zl = Yl;
% Loc.Xfs = Xfs;
% Loc.Zfs = Yfs;
% Loc.Xrs = Xrs;
% Loc.Zrs = Yrs;
Loc = nan;


Failure.U.tension =  g_tension_u;
Failure.U.compression = g_compression_u;
Failure.U.buckling = g_buckling_u;

Failure.L.tension =  g_tension_l;
Failure.L.compression = g_compression_l;
Failure.L.buckling = g_buckling_l;

Failure.FS.shear =  g_shear_fs;
Failure.FS.buckling = g_buckling_fs;

Failure.RS.shear =  g_shear_rs;
Failure.RS.buckling = g_buckling_rs;

%% Plot

if plotoption ==1
    figure
    hold on
    axis equal
    
    for i=1:size(Loc.Y,1)
        plot3(Loc.Xu(i,:),Loc.Y(i)*ones(1,10),Loc.Zu(i,:));
        plot3(Loc.Xl(i,:),Loc.Y(i)*ones(1,10),Loc.Zl(i,:));
        plot3(Loc.Xfs(i,:),Loc.Y(i)*ones(1,5),Loc.Zfs(i,:));
        plot3(Loc.Xrs(i,:),Loc.Y(i)*ones(1,5),Loc.Zrs(i,:));
    end
end

end

