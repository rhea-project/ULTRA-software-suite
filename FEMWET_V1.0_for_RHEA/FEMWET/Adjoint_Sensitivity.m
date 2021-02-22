function [Dg_Dx, dU_dX, dMa_dX, dCD_dX] = Adjoint_Sensitivity(AC,X,D,Z,dg_dx,dK_dx,dF_dx,dAIC_dx,dV_dx,dL_dx,U,Gamma,dM1_dx, dM2_dx, dM_dGamma,D2,dK2_dx,dF2_dx,dAIC2_dx,dV2_dx,dL2_dx,U2,Gamma2,dCL_dx, dCD_dX, Zcd, Drag, Aileron_Eff,dWdes_dx) %dA_da,dAIC_da, d2AIC_dxda, d2V_dxda,si_m)


nElm     = AC.FEM.nElm;
nSkin    = AC.FEM.nSkin;
nSpar    = AC.FEM.nSpar;
ActiveDOFs = AC.FEM.ActiveDOFs;
nChord   = AC.Q3D.Grid.nChord;
nDOF = AC.FEM.nDOF;

nX = length(X);
nZ = nElm*nChord+nDOF-6+1;
n2d = AC.Q3D.N2D;

%% dS_dx, dA_dx, dW_dx

dS_dx = zeros(nDOF-6,nX);
dA_dx = zeros(nElm*nChord,nX);
dW_dx = zeros(1,nX);
% dCL_dx = zeros(AC.Grid.N2D,nX);

for k = 1:nX
    df_dx(1:nDOF-6,1) = dF_dx(ActiveDOFs,k);
    dk_dx(1:nDOF-6,1:nDOF-6) = dK_dx(ActiveDOFs,ActiveDOFs,k);
   
    da_dx(1:nElm*nChord,1:nElm*nChord) = dAIC_dx(:,:,k);
    dv_dx(1:nElm*nChord,1) = dV_dx(:,k); 

    dS_dx(:,k) = dk_dx*U(ActiveDOFs) - df_dx; 
    dA_dx(:,k) = da_dx*Gamma - dv_dx;
    dW_dx(:,k) = dL_dx(k) - dWdes_dx(k);
%     dCL_dx(:,k) = dCL_dx(:,k);
end

%% dU_dX
if Drag ==1
    Zu = [zeros(nElm*nChord,nDOF-6); eye(nDOF-6,nDOF-6); zeros(1,nDOF-6); zeros(n2d,nDOF-6);];
    
    Adj_u = D'\Zu;
    si_u =  Adj_u(1:nElm*nChord,:);   phi_u = Adj_u(nElm*nChord+1:nElm*nChord+nDOF-6,:);   lambda_u = Adj_u(nElm*nChord+nDOF-6+1,:);  
    psi_u = Adj_u(nElm*nChord+nDOF-6+2:end,:);

    dU_dX = - si_u'*dA_dx  - phi_u'*dS_dx - lambda_u'*dW_dx - psi_u'*dCL_dx;
    dU_dX = [zeros(6,nX); dU_dX];

else
    Zu = [zeros(nElm*nChord,nDOF-6); eye(nDOF-6,nDOF-6); zeros(1,nDOF-6)];
    Adj_u = D'\Zu;
    si_u =  Adj_u(1:nElm*nChord,:);   phi_u = Adj_u(nElm*nChord+1:nElm*nChord+nDOF-6,:);   lambda_u = Adj_u(nElm*nChord+nDOF-6+1,:);  
    dU_dX = - si_u'*dA_dx  - phi_u'*dS_dx - lambda_u'*dW_dx;
    dU_dX = [zeros(6,nX); dU_dX];
end




%% upper and lower stress

for i=1:nElm
    for j=1:nSkin

        pg_px_tu(1:nX,1) = dg_dx.U.tension(i,j,:);
        pg_px_cu(1:nX,1) = dg_dx.U.compression(i,j,:);
        pg_px_bu(1:nX,1) = dg_dx.U.buckling(i,j,:);
        
        pg_px_tl(1:nX,1) = dg_dx.L.tension(i,j,:);
        pg_px_cl(1:nX,1) = dg_dx.L.compression(i,j,:);
        pg_px_bl(1:nX,1) = dg_dx.L.buckling(i,j,:);
        
        z_tu(1:nZ,1) = Z.U.tension(i,:,j);
        z_cu(1:nZ,1) = Z.U.compression(i,:,j);
        z_bu(1:nZ,1) = Z.U.buckling(i,:,j);
        
        z_tl(1:nZ,1) = Z.L.tension(i,:,j);
        z_cl(1:nZ,1) = Z.L.compression(i,:,j);
        z_bl(1:nZ,1) = Z.L.buckling(i,:,j);


        if Drag ==1
            Adj_cu = D'\[z_cu; zeros(n2d,1)];
            Adj_bu = D'\[z_bu; zeros(n2d,1)];
            Adj_tu = D'\[z_tu; zeros(n2d,1)];
            
            Adj_cl = D'\[z_cl; zeros(n2d,1)];
            Adj_bl = D'\[z_bl; zeros(n2d,1)];
            Adj_tl = D'\[z_tl; zeros(n2d,1)];

            si_tu = Adj_tu(1:nElm*nChord,1);   phi_tu = Adj_tu(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_tu = Adj_tu(nElm*nChord+nDOF-6+1,1);   psi_tu = Adj_tu(nElm*nChord+nDOF-6+2:end);
            si_cu = Adj_cu(1:nElm*nChord,1);   phi_cu = Adj_cu(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_cu = Adj_cu(nElm*nChord+nDOF-6+1,1);   psi_cu = Adj_cu(nElm*nChord+nDOF-6+2:end);
            si_bu = Adj_bu(1:nElm*nChord,1);   phi_bu = Adj_bu(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_bu = Adj_bu(nElm*nChord+nDOF-6+1,1);   psi_bu = Adj_bu(nElm*nChord+nDOF-6+2:end);

            si_cl = Adj_cl(1:nElm*nChord,1);   phi_cl = Adj_cl(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_cl = Adj_cu(nElm*nChord+nDOF-6+1,1);   psi_cl = Adj_cl(nElm*nChord+nDOF-6+2:end);
            si_bl = Adj_bl(1:nElm*nChord,1);   phi_bl = Adj_bl(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_bl = Adj_bu(nElm*nChord+nDOF-6+1,1);   psi_bl = Adj_bl(nElm*nChord+nDOF-6+2:end);
            si_tl = Adj_tl(1:nElm*nChord,1);   phi_tl = Adj_tl(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_tl = Adj_tl(nElm*nChord+nDOF-6+1,1);   psi_tl = Adj_tl(nElm*nChord+nDOF-6+2:end);


            for k = 1:nX
                Dg_Dx.U.tension(i,j,k)     = pg_px_tu(k) - si_tu'*dA_dx(:,k) - phi_tu'*dS_dx(:,k) - lambda_tu'*dW_dx(:,k) - psi_tu'*dCL_dx(:,k) ;
                Dg_Dx.U.compression(i,j,k) = pg_px_cu(k) - si_cu'*dA_dx(:,k) - phi_cu'*dS_dx(:,k) - lambda_cu'*dW_dx(:,k) - psi_cu'*dCL_dx(:,k) ;
                Dg_Dx.U.buckling(i,j,k)    = pg_px_bu(k) - si_bu'*dA_dx(:,k) - phi_bu'*dS_dx(:,k) - lambda_bu'*dW_dx(:,k) - psi_bu'*dCL_dx(:,k) ; 
                  
                Dg_Dx.L.tension(i,j,k)     = pg_px_tl(k) - si_tl'*dA_dx(:,k) - phi_tl'*dS_dx(:,k) - lambda_tl'*dW_dx(:,k) - psi_tl'*dCL_dx(:,k) ; 
                Dg_Dx.L.compression(i,j,k) = pg_px_cl(k) - si_cl'*dA_dx(:,k) - phi_cl'*dS_dx(:,k) - lambda_cl'*dW_dx(:,k) - psi_cl'*dCL_dx(:,k) ;
                Dg_Dx.L.buckling(i,j,k)    = pg_px_bl(k) - si_bl'*dA_dx(:,k) - phi_bl'*dS_dx(:,k) - lambda_bl'*dW_dx(:,k) - psi_bl'*dCL_dx(:,k) ;
            end
        else
            
            Adj_tu = D'\z_tu;
            Adj_cu = D'\z_cu;
            Adj_bu = D'\z_bu;
        
            Adj_tl = D'\z_tl;
            Adj_cl = D'\z_cl;
            Adj_bl = D'\z_bl;
 
            si_tu = Adj_tu(1:nElm*nChord,1);   phi_tu = Adj_tu(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_tu = Adj_tu(nElm*nChord+nDOF-6+1,1);
            si_cu = Adj_cu(1:nElm*nChord,1);   phi_cu = Adj_cu(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_cu = Adj_cu(nElm*nChord+nDOF-6+1,1);
            si_bu = Adj_bu(1:nElm*nChord,1);   phi_bu = Adj_bu(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_bu = Adj_bu(nElm*nChord+nDOF-6+1,1);

            si_tl = Adj_tl(1:nElm*nChord,1);   phi_tl = Adj_tl(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_tl = Adj_tl(nElm*nChord+nDOF-6+1,1);
            si_cl = Adj_cl(1:nElm*nChord,1);   phi_cl = Adj_cl(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_cl = Adj_cl(nElm*nChord+nDOF-6+1,1);
            si_bl = Adj_bl(1:nElm*nChord,1);   phi_bl = Adj_bl(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_bl = Adj_bl(nElm*nChord+nDOF-6+1,1);

            for k = 1:nX
                Dg_Dx.U.tension(i,j,k)     = pg_px_tu(k) - si_tu'*dA_dx(:,k) - phi_tu'*dS_dx(:,k) - lambda_tu'*dW_dx(:,k);
                Dg_Dx.U.compression(i,j,k) = pg_px_cu(k) - si_cu'*dA_dx(:,k) - phi_cu'*dS_dx(:,k) - lambda_cu'*dW_dx(:,k);
                Dg_Dx.U.buckling(i,j,k)    = pg_px_bu(k) - si_bu'*dA_dx(:,k) - phi_bu'*dS_dx(:,k) - lambda_bu'*dW_dx(:,k);

                Dg_Dx.L.tension(i,j,k)     = pg_px_tl(k) - si_tl'*dA_dx(:,k) - phi_tl'*dS_dx(:,k) - lambda_tl'*dW_dx(:,k);
                Dg_Dx.L.compression(i,j,k) = pg_px_cl(k) - si_cl'*dA_dx(:,k) - phi_cl'*dS_dx(:,k) - lambda_cl'*dW_dx(:,k);
                Dg_Dx.L.buckling(i,j,k)    = pg_px_bl(k) - si_bl'*dA_dx(:,k) - phi_bl'*dS_dx(:,k) - lambda_bl'*dW_dx(:,k);    
            end
        end
    end
end
       

%% Spars stress

for i=1:nElm
    for j=1:nSpar

        pg_px_sfs(1:nX,1) = dg_dx.FS.shear(i,j,:);
        pg_px_bfs(1:nX,1) = dg_dx.FS.buckling(i,j,:);

        pg_px_srs(1:nX,1) = dg_dx.RS.shear(i,j,:);
        pg_px_brs(1:nX,1) = dg_dx.RS.buckling(i,j,:);


        z_sfs(1:nZ,1) = Z.FS.shear(i,:,j);
        z_bfs(1:nZ,1) = Z.FS.buckling(i,:,j);

        z_srs(1:nZ,1) = Z.RS.shear(i,:,j);
        z_brs(1:nZ,1) = Z.RS.buckling(i,:,j);

        
        if Drag ==1
        
            Adj_sfs = D'\[z_sfs; zeros(n2d,1)];
            Adj_bfs = D'\[z_bfs; zeros(n2d,1)];

            Adj_srs = D'\[z_srs; zeros(n2d,1)];
            Adj_brs = D'\[z_brs; zeros(n2d,1)];

            si_sfs = Adj_sfs(1:nElm*nChord,1);   phi_sfs = Adj_sfs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_sfs = Adj_sfs(nElm*nChord+nDOF-6+1,1);    psi_sfs = Adj_sfs(nElm*nChord+nDOF-6+2:end);
            si_bfs = Adj_bfs(1:nElm*nChord,1);   phi_bfs = Adj_bfs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_bfs = Adj_bfs(nElm*nChord+nDOF-6+1,1);    psi_bfs = Adj_bfs(nElm*nChord+nDOF-6+2:end);

            si_srs = Adj_srs(1:nElm*nChord,1);   phi_srs = Adj_srs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_srs = Adj_srs(nElm*nChord+nDOF-6+1,1);    psi_srs = Adj_srs(nElm*nChord+nDOF-6+2:end);
            si_brs = Adj_brs(1:nElm*nChord,1);   phi_brs = Adj_brs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_brs = Adj_brs(nElm*nChord+nDOF-6+1,1);    psi_brs = Adj_brs(nElm*nChord+nDOF-6+2:end);

            for k = 1:nX
                Dg_Dx.FS.shear(i,j,k)     = pg_px_sfs(k) - si_sfs'*dA_dx(:,k) - phi_sfs'*dS_dx(:,k) - lambda_sfs'*dW_dx(:,k) - psi_sfs'*dCL_dx(:,k);
                Dg_Dx.FS.buckling(i,j,k)  = pg_px_bfs(k) - si_bfs'*dA_dx(:,k) - phi_bfs'*dS_dx(:,k) - lambda_bfs'*dW_dx(:,k) - psi_bfs'*dCL_dx(:,k);

                Dg_Dx.RS.shear(i,j,k)     = pg_px_srs(k) - si_srs'*dA_dx(:,k) - phi_srs'*dS_dx(:,k) - lambda_srs'*dW_dx(:,k) - psi_srs'*dCL_dx(:,k);  
                Dg_Dx.RS.buckling(i,j,k)  = pg_px_brs(k) - si_brs'*dA_dx(:,k) - phi_brs'*dS_dx(:,k) - lambda_brs'*dW_dx(:,k) - psi_brs'*dCL_dx(:,k);
            end
        else
            Adj_sfs = D'\z_sfs;
            Adj_bfs = D'\z_bfs;

            Adj_srs = D'\z_srs;
            Adj_brs = D'\z_brs;

            si_sfs = Adj_sfs(1:nElm*nChord,1);   phi_sfs = Adj_sfs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_sfs = Adj_sfs(nElm*nChord+nDOF-6+1,1);    
            si_bfs = Adj_bfs(1:nElm*nChord,1);   phi_bfs = Adj_bfs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_bfs = Adj_bfs(nElm*nChord+nDOF-6+1,1);   

            si_srs = Adj_srs(1:nElm*nChord,1);   phi_srs = Adj_srs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_srs = Adj_srs(nElm*nChord+nDOF-6+1,1);    
            si_brs = Adj_brs(1:nElm*nChord,1);   phi_brs = Adj_brs(nElm*nChord+1:nElm*nChord+nDOF-6);   lambda_brs = Adj_brs(nElm*nChord+nDOF-6+1,1);  

            for k = 1:nX
                Dg_Dx.FS.shear(i,j,k)     = pg_px_sfs(k) - si_sfs'*dA_dx(:,k) - phi_sfs'*dS_dx(:,k) - lambda_sfs'*dW_dx(:,k);
                Dg_Dx.FS.buckling(i,j,k)  = pg_px_bfs(k) - si_bfs'*dA_dx(:,k) - phi_bfs'*dS_dx(:,k) - lambda_bfs'*dW_dx(:,k);

                Dg_Dx.RS.shear(i,j,k)     = pg_px_srs(k) - si_srs'*dA_dx(:,k) - phi_srs'*dS_dx(:,k) - lambda_srs'*dW_dx(:,k); 
                Dg_Dx.RS.buckling(i,j,k)  = pg_px_brs(k) - si_brs'*dA_dx(:,k) - phi_brs'*dS_dx(:,k) - lambda_brs'*dW_dx(:,k);
            end 
            
        end
    end
end


%% **************** Sensitivity of Aileron efficiency:  deta_dx *******************

if Aileron_Eff ==1
    dS2_dx = zeros(nDOF-6,nX);
    dA2_dx = zeros(nElm*nChord,nX);
    dW2_dx = zeros(1,nX);

    for k = 1:nX
        df_dx(1:nDOF-6,1) = dF2_dx(ActiveDOFs,k);
        dk_dx(1:nDOF-6,1:nDOF-6) = dK2_dx(ActiveDOFs,ActiveDOFs,k);

        da_dx(1:nElm*nChord,1:nElm*nChord) = dAIC2_dx(:,:,k);
        dv_dx(1:nElm*nChord,1) = dV2_dx(:,k); 

        dS2_dx(:,k) = dk_dx*U2(ActiveDOFs) - df_dx; 
        dA2_dx(:,k) = da_dx*Gamma2 - dv_dx;
        dW2_dx(:,k) = dL2_dx(k) - dWdes_dx(k);
    end


    Zm = [dM_dGamma; zeros(length(U)-6,1); 0]; % zeros(n2d,1)];

    Adj_m1 = D'\Zm;
    si_m1 =  Adj_m1(1:nElm*nChord,:);   phi_m1 = Adj_m1(nElm*nChord+1:nElm*nChord+nDOF-6,:);   lambda_m1 = Adj_m1(nElm*nChord+nDOF-6+1,:); % psi_m1 = Adj_m1(nElm*nChord+nDOF-6+2:end,:);
    dM1_dX = dM1_dx - si_m1'*dA_dx - phi_m1'*dS_dx - lambda_m1'*dW_dx; % - psi_m1'*dCL_dx;

    % Zm2 = [dM_dGamma; zeros(length(U)-6,1); 0];
    Adj_m2 = D2'\Zm;
    si_m2 =  Adj_m2(1:nElm*nChord,:);   phi_m2 = Adj_m2(nElm*nChord+1:nElm*nChord+nDOF-6,:);   lambda_m2 = Adj_m2(nElm*nChord+nDOF-6+1,:); % psi_m2 = Adj_m2(nElm*nChord+nDOF-6+2:end,:);
    dM2_dX = dM2_dx - si_m2'*dA2_dx - phi_m2'*dS2_dx - lambda_m2'*dW2_dx; % - psi_m2'*dCL_dx; 

    dMa_dX = (dM2_dX-dM1_dX)/(1e-6);

else
    dMa_dX = nan;
end


%% dCD_dX
if Drag ==1
    Adj_cd = D'\Zcd;
    si_cd =  Adj_cd(1:nElm*nChord,:);   phi_cd = Adj_cd(nElm*nChord+1:nElm*nChord+nDOF-6,:);   lambda_cd = Adj_cd(nElm*nChord+nDOF-6+1,:);  psi_cd = Adj_cd(nElm*nChord+nDOF-6+2:end,:);

    dCD_dX = dCD_dX - si_cd'*dA_dx - phi_cd'*dS_dx - lambda_cd'*dW_dx - psi_cd'*dCL_dx;
else
    dCD_dX = nan;
end



