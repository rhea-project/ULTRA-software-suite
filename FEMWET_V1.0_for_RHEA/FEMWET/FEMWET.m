function [Wwing, Failure, U, eta_a, CD, WL, Vf, dWw_dX, Dg_Dx,dU_dX, dMa_dX, dCD_dX, dWL_dX, dVf_dX, CDi, CDp, CDf, dM_da,dM_da0]...
    =FEMWET(AC,T,G,P,Case,Adjoint_Switch,Drag,Aileron_Eff,DV,Wf,MTOW)

%% Adjoint 

global Adjoint SF Restart

SF = 1.5;
Restart = 1;


if Drag==1
    Aileron_Eff =0;
end

%% Design variables

if DV == 1
    X  = T;
    DP = [G; P];
elseif DV == 2
    X  = G;
    DP = [T; P];
elseif DV == 3
    X  = P;
    DP = [T; G];
elseif DV ==4
    X  = [G; P];
    DP = T;
elseif DV ==5
    X  = [T; G];
    DP = P;
elseif DV ==6
    X  = [T; P];
    DP = G;
elseif DV ==7
    X  = [T; G; P];
    DP = [];
end

% adding fuel weight and MTOW as surrogate variables to X
X = [X; Wf; MTOW]; 

%% Load Cases

H = AC.Load{Case}.H;
Mach = AC.Load{Case}.Mach;
n = AC.Load{Case}.n;

[T, a, ~, rho] = atmosisa(H);
Vinf = Mach*a; 

    
if Case ==5
    AC.Aileron.Geom(end) = 0.01;  % for aileron effectiveness calculation
end


AC.Aero.Mach = Mach;
AC.Aero.Vinf = Vinf;        % speed
AC.Aero.rho = rho;                    % air density
AC.Aero.mu  = 1.716e-5*(273.15+110.3)/(T+110.3)*(T/273.15)^(3/2);   % Sutherland's law

%% Weights

dWdes_dx = zeros(size(X));
dWfload_dx = zeros(size(X));

if Case == 1 || Case == 2 || Case == 3 % || Case == 4 % for Wdes = MTOW
    Wdes = n*MTOW*9.81;
    dWdes_dx(end) = n*9.81;
    
    Wfload = Wf*9.81;
    dWfload_dx(end-1) = 9.81;

elseif Case ==4             % zero fuel weight for gust
    Wdes = n*(MTOW-Wf)*9.81;
    dWdes_dx(end) =n*9.81;
    dWdes_dx(end-1) = -n*9.81; 
    
    Wfload = 0;
    
elseif Case == 5 || Case == 6 % for Wdes = cruise weight
    Wdes = sqrt(MTOW*(MTOW-Wf))*9.81;
    dWdes_dx(end) = (2*MTOW-Wf)*0.5*(MTOW^2-MTOW*Wf)^(-0.5)*9.81;
    dWdes_dx(end-1) = -MTOW*0.5*(MTOW^2-MTOW*Wf)^(-0.5)*9.81;
    
    Wfload = MTOW*9.81 - Wdes;
    dWfload_dx(end) = 9.81- (2*MTOW-Wf)*0.5*(MTOW^2-MTOW*Wf)^(-0.5)*9.81;
    dWfload_dx(end-1) = MTOW*0.5*(MTOW^2-MTOW*Wf)^(-0.5)*9.81;
end
    
   

%% Connectivity Matrix

e = 1:AC.FEM.nElm;
AC.FEM.Conc(:,1) = e;
AC.FEM.Conc(:,2) = e+1;

%% Fixed DOFs

AC.FEM.FixedDOFs = 1:6; 
AC.FEM.ActiveDOFs=setdiff(1:AC.FEM.nDOF,AC.FEM.FixedDOFs);


%% Nodes

if Adjoint_Switch ==0
    [cNode, Xea, Yea] = Node_Generation(AC,X,DV,DP);
    dcNode = zeros(size(cNode,1),size(X,1));
elseif Adjoint_Switch ==1
    [dcNode, dXea, dYea] =Node_Generation(AC,gradientinit(X),DV,DP);

    cNode = dcNode.x;
    dcNode = dcNode.dx;

    Xea = dXea.x;
    dXea = dXea.dx;

    Yea = dYea.x;
    dYea = dYea.dx;
end


%% Fuel volume

dVf = fuel_volume(gradientinit(X),AC,DV,DP);
Vf = dVf.x;
dVf_dX = dVf.dx;


%% Newton method

Adjoint = 0;

Gamma = zeros(AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord,1);
U=zeros(AC.FEM.nDOF,1);
Alpha = 1*pi/180;
Alpha_i = zeros(AC.Q3D.N2D,1);    

nS = sum(AC.Q3D.Grid.DnSpan);
    
for iter = 1:20
    % derivatives
    % structural part
    [K, F, ~, ~, ~, ~, dF_dL] = Matrices(AC,AC.Aero.rho*AC.Aero.Vinf*Gamma/sqrt(1-AC.Aero.Mach^2),X,cNode,Xea,Yea,0,0,0,n,Vf,dVf_dX,Wfload,dWfload_dx,DV,DP);
    dF_dGamma = dF_dL*(AC.Aero.rho*AC.Aero.Vinf/sqrt(1-AC.Aero.Mach^2));

    [dA, dV, dV_dA, dB] = AIC_RHS_Aileron(X,AC,Alpha,gradientinit(U),Xea,DV,DP);
    
    A = dA.x;
    V = dV.x;
    dA_dU = dA.dx;
    dV_dU = dV.dx;
    dV_dAlpha = dV_dA.x;
    
    B = dB.x;
    dB_dU = dB.dx;
    dB_dAlpha = 0;
%     Wind = B*Gamma;
    G = Gamma(end-nS+2:end);
    G = [G(1:end-1)-G(2:end); G(end)]/sqrt(1-AC.Aero.Mach^2);
    Wind = B*G;
    
    dA_dAlpha = 0;
    dGamma_dAlpha = A\(dV_dAlpha-dA_dAlpha*Gamma);
    %     dWind_dAlpha  = dB_dAlpha*Gamma + B*dGamma_dAlpha;
    dG_dAlpha = dGamma_dAlpha(end-nS+2:end);
    dG_dAlpha = [dG_dAlpha(1:end-1)-dG_dAlpha(2:end); dG_dAlpha(end)]/sqrt(1-AC.Aero.Mach^2);
    dWind_dAlpha  = dB_dAlpha*G + B*dG_dAlpha;
      
    
    dAV_dU = zeros(AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord,AC.FEM.nDOF);
    dGamma_dU = zeros(AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord,AC.FEM.nDOF);
    dWind_dU = zeros(nS-1,AC.FEM.nDOF);
    
    
    for i=1:length(U)
        dAV_dU(:,i) = dA_dU(:,:,i)*Gamma - dV_dU(:,i);
        dGamma_dU(:,i) = A\(dV_dU(:,i)-dA_dU(:,:,i)*Gamma);
        dG_dU = [dGamma_dU(end-nS+2:end-1,i)-dGamma_dU(end-nS+3:end,i); dGamma_dU(end,i)]/sqrt(1-AC.Aero.Mach^2);
        dWind_dU(:,i) = dB_dU(:,:,i)*G + B*dG_dU;
    end


    [Ltot, ~, ~, ~, clvlm, ~, ~, ~, dCDivlm, dCDivlm_dWind] = Lift_VLM(X,AC,gradientinit(Gamma),Wind,DV,DP);
    
    Lift = Ltot.x;
    dL_dGamma = full(Ltot.dx);    
    Clvlm = clvlm.x;
    dClvlm_dGamma = clvlm.dx;
    
    CDivlm = dCDivlm.x;
    dCDivlm_dWind = dCDivlm_dWind.x;
    dCDivlm_dGamma = dCDivlm.dx; 
        
    dCDivlm_dAlpha = dCDivlm_dGamma*dGamma_dAlpha + dCDivlm_dWind*dWind_dAlpha;   
    dCDivlm_dU     = dCDivlm_dGamma*dGamma_dU + dCDivlm_dWind*dWind_dU;      


    
    if Drag ==1
        [Cl2d, CD, dCl2d_dX, pCD_pX, dCl2d_dAlpha, dCD_dAlpha, dCl2d_dAlpha_i, dCD_dAlpha_i, dCl2d_dU, dCD_dU, ~, CDp, CDf] =Aero_coef(AC,X,Alpha,Alpha_i,U,iter,1,cNode,dcNode,DV,DP);  

        dCL_dGamma = -dClvlm_dGamma;
        dCL_dU     = dCl2d_dU;
        dCL_dAlpha = dCl2d_dAlpha;
        dCL_dAlpha_i = dCl2d_dAlpha_i;

     % for CDi from VLM
        CDi = CDivlm;
        CD = CD + CDivlm;
        dCD_dAlpha = dCD_dAlpha + dCDivlm_dAlpha;
        % dCDivlm_dAlpha_i = 0;
        
        dCD_dU = dCD_dU + dCDivlm_dU;
       % ******
        
        
        
        D = [A                                   dAV_dU(:,AC.FEM.ActiveDOFs)               -dV_dAlpha                  zeros(length(dV_dAlpha),length(Alpha_i));
            -dF_dGamma(AC.FEM.ActiveDOFs,:)      K(AC.FEM.ActiveDOFs,AC.FEM.ActiveDOFs)    zeros(AC.FEM.nDOF-6,1)      zeros(AC.FEM.nDOF-6,length(Alpha_i));
             dL_dGamma                           zeros(1,AC.FEM.nDOF-6)                    zeros(1,1)                  zeros(1,length(Alpha_i));
             dCL_dGamma                          dCL_dU(:,AC.FEM.ActiveDOFs)               dCL_dAlpha                  full(dCL_dAlpha_i)];



        RHS = -[A*Gamma-V;  K(AC.FEM.ActiveDOFs,AC.FEM.ActiveDOFs)*U(AC.FEM.ActiveDOFs)-F(AC.FEM.ActiveDOFs); Lift-Wdes; Cl2d-Clvlm];

        deltas = D\RHS;

        dGamma = deltas(1:AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord, 1);
        du = deltas(AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+1:AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+AC.FEM.nDOF-6, 1);
        dalpha = deltas(AC.FEM.nDOF-6+AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+1, 1);
        dalpha_i = deltas(AC.FEM.nDOF-6+AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+2:end,1);

        U(AC.FEM.ActiveDOFs) = U(AC.FEM.ActiveDOFs) + du;
        Gamma = Gamma + dGamma;
        Alpha = Alpha + dalpha;
        Alpha_i = Alpha_i + dalpha_i;
    
    else
        D = [A                               dAV_dU(:,AC.FEM.ActiveDOFs)               -dV_dAlpha                 
        -dF_dGamma(AC.FEM.ActiveDOFs,:)      K(AC.FEM.ActiveDOFs,AC.FEM.ActiveDOFs)    zeros(AC.FEM.nDOF-6,1)     
         dL_dGamma                           zeros(1,AC.FEM.nDOF-6)                    zeros(1,1)               ];
      
 
        RHS = -[A*Gamma-V;  K(AC.FEM.ActiveDOFs,AC.FEM.ActiveDOFs)*U(AC.FEM.ActiveDOFs)-F(AC.FEM.ActiveDOFs); Lift-Wdes];

        deltas = D\RHS;

        dGamma = deltas(1:AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord, 1);
        du = deltas(AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+1:AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+AC.FEM.nDOF-6, 1);
        dalpha = deltas(AC.FEM.nDOF-6+AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+1, 1);


        U(AC.FEM.ActiveDOFs) = U(AC.FEM.ActiveDOFs) + du;
        Gamma = Gamma + dGamma;
        Alpha = Alpha + dalpha;

        CD = nan; CDi = nan;  CDp = nan;  CDf = nan;
    end

    Error(iter) = norm(deltas, inf);

    if Error(iter) < 1E-4 || isnan(Error(iter))
        break;
    elseif Error(iter) > 1E10 
        fprintf (1, 'Wing is diverging!');
        error('Wing is diverging!');
    end

end



if isnan(Error(iter))
   Failure= nan;
   U = nan;
   eta_a = nan;
   CD = nan;
   Dg_Dx = nan;
   dU_dX = nan;
   dMa_dX = nan;
   dCD_dX = nan;
   CDi = nan;
   CDp = nan;
   CDf = nan;
   return
end

if Error(iter) > 1e-4
%     error('The Newton iteration is not converged, increase the maximum number of iterations');
    fprintf (1,'Newton iteration is not converged, increase the maximum number of iterations');
    Failure = nan;
    U = nan;
    eta_a = nan;
    CD = nan;
    Dg_Dx = nan;
    dU_dX = nan;
    dMa_dX = nan;
    dCD_dX = nan;
    CDi = nan;
    CDp = nan;
    CDf = nan;
    dM_da = nan;
    dM_da0 = nan;
    return
end

%% Finding M2

if Aileron_Eff ==1
    
    AC.Aileron.Geom(end) = AC.Aileron.Geom(end) + 1e-6*180/pi;

    Gamma2 = Gamma;
    U2 = U;
    Alpha2 = Alpha;
    % Alpha_i2 = Alpha_i;

    clear Error

    for iter = 1:10
        % derivatives
        [K, F, ~, ~, ~, ~, dF_dL] = Matrices(AC,AC.Aero.rho*AC.Aero.Vinf*Gamma2/sqrt(1-AC.Aero.Mach^2),X,cNode,Xea,Yea,0,0,0,n,Vf,dVf_dX,Wfload,dWfload_dx,DV,DP);
        dF_dGamma =  dF_dL*(AC.Aero.rho*AC.Aero.Vinf/sqrt(1-AC.Aero.Mach^2));

        [dA, dV, dV_dA] = AIC_RHS_Aileron(X,AC,Alpha2,gradientinit(U2),Xea,DV,DP);

        A = dA.x;
        V = dV.x;
        dA_dU = dA.dx;
        dV_dU = dV.dx;
        dV_dAlpha = dV_dA.x;
        dAV_dU = zeros(AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord,AC.FEM.nDOF);
        for i=1:length(U)
            dAV_dU(:,i) = dA_dU(:,:,i)*Gamma2 - dV_dU(:,i);
        end

%         B = dB.x;
%         Wind2 = B*Gamma2;
        Wind2 = 0;

        Ltot = Lift_VLM(X,AC,gradientinit(Gamma2),Wind2,DV,DP);
        Lift = Ltot.x;
        dL_dGamma = full(Ltot.dx);    


        D2 = [A                                  dAV_dU(:,AC.FEM.ActiveDOFs)              -dV_dAlpha                 
            -dF_dGamma(AC.FEM.ActiveDOFs,:)      K(AC.FEM.ActiveDOFs,AC.FEM.ActiveDOFs)    zeros(AC.FEM.nDOF-6,1)      
             dL_dGamma                           zeros(1,AC.FEM.nDOF-6)                    zeros(1,1)               ];   


        RHS = -[A*Gamma2-V;  K(AC.FEM.ActiveDOFs,AC.FEM.ActiveDOFs)*U2(AC.FEM.ActiveDOFs)-F(AC.FEM.ActiveDOFs); Lift-Wdes];

        deltas = D2\RHS;

        dGamma = deltas(1:AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord, 1);
        du = deltas(AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+1:AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+AC.FEM.nDOF-6, 1);
        dalpha = deltas(AC.FEM.nDOF-6+AC.Q3D.Grid.nSpan*AC.Q3D.Grid.nChord+1, 1);
  
        U2(AC.FEM.ActiveDOFs) = U2(AC.FEM.ActiveDOFs) + du;
        Gamma2 = Gamma2 + dGamma;
        Alpha2 = Alpha2 + dalpha;
   

        Error(iter) = norm(deltas, inf); 
        if Error(iter) < 1E-6
            break;
        elseif Error(iter) > 1E10 || isnan(Error(iter))
            fprintf (1, 'Wing is diverging!');
            stop
        end

    end

    if Error(iter) > 1e-6
        error('The Newton iteration is not converged, increase the maximum number of iterations');
    end

    AC.Aileron.Geom(end) = AC.Aileron.Geom(end) - 1e-6*180/pi;

    % ************** Aileron efficiency **************

    [eta_a, dM_da, dM_da0] = Aileron_Efficiency(X,AC,D,Wdes,Alpha,U,Xea,AC.Aileron.Geom(end)*pi/180,DV,DP);

else
   eta_a = nan; 
   dM_da = nan;
   dM_da0 = nan;
end

%% Stress

Failure = Stress(AC,U,X,cNode,Xea,Yea,0,0,0,DV,DP);

%% Weight

Wwing = Weight(AC,X,DV,DP);
WL = WingLoading(X,AC,DV,DP);

%% Fuel tank volume

Vf = fuel_volume(X,AC,DV,DP);

%% plotting
% 
% u = U(1:6:end);    % displacement in x-direction
% v = U(2:6:end);    % displacement in y-direction 
% w = U(3:6:end);    % displacement in z-direction (bending)
% 
% px = U(4:6:end)*180/pi;   % deflection angle x-axis in degree
% py = U(5:6:end)*180/pi;   % deflection angle y-axis (twist)
% pz = U(6:6:end)*180/pi;   % deflection angle z-axis 
% 
% Wing_plot2(AC,AC.FEM.nNode,cNode,u,v,w,px,py,pz,0);


%% ******************************* Adjoint analysis ***********************************


if Adjoint_Switch ==1

Adjoint = 1;


%% Aerodynamic

Load = AC.Aero.rho*AC.Aero.Vinf*Gamma/sqrt(1-AC.Aero.Mach^2) ;
[dA, dRHS] = AIC_RHS_Aileron(gradientinit(X),AC,Alpha,U,Xea,DV,DP);
dA_dx = dA.dx;
dRHS_dx = dRHS.dx;

[dA, dRHS] = AIC_RHS_Aileron(X,AC,Alpha,U,gradientinit(Xea),DV,DP);
dA_dXea = dA.dx;
dRHS_dXea = dRHS.dx;

for i=1:length(X)
    dA_dX(:,:,i) = dA_dx(:,:,i) + dA_dXea*dXea(i);
    dRHS_dX(:,i) = dRHS_dx(:,i) + dRHS_dXea*dXea(i);
end       

[dLtot, dM,~,~,dClvlm] = Lift_VLM(gradientinit(X),AC,Gamma,Wind,DV,DP);
dLtot_dX = dLtot.dx;
dM1_dX  = dM.dx;
dClvlm_dX = dClvlm.dx;

% **** for dM_da ****
if Aileron_Eff ==1
    Load2 =  AC.Aero.rho*AC.Aero.Vinf*Gamma2/sqrt(1-AC.Aero.Mach^2) ;
    [dA2, dRHS2] = AIC_RHS_Aileron(gradientinit(X),AC,Alpha2,U2,Xea,DV,DP,AC.Aileron.Geom(end)*pi/180+1e-6);
    dA2_dx = dA2.dx;
    dRHS2_dx = dRHS2.dx;

    [dA2, dRHS2] =  AIC_RHS_Aileron(X,AC,Alpha2,U2,gradientinit(Xea),DV,DP,AC.Aileron.Geom(end)*pi/180+1e-6);
    dA2_dXea = dA2.dx;
    dRHS2_dXea = dRHS2.dx;

    for i=1:length(X)
        dA2_dX(:,:,i) = dA2_dx(:,:,i) + dA2_dXea*dXea(i);
        dRHS2_dX(:,i) = dRHS2_dx(:,i) + dRHS2_dXea*dXea(i);
    end       

    [dLtot2, dM2, ~, ~, ~, ~,y, dy] =  Lift_VLM(gradientinit(X),AC,Gamma2,Wind2,DV,DP);
    dLtot2_dX = dLtot2.dx;
    dM2_dX  = dM2.dx;
    y = y.x;
    dy = dy.x;

    dM_dGamma = zeros(size(Gamma));
    dM_dGamma(end-AC.Q3D.Grid.nSpan+1:end) =  AC.Aero.rho*AC.Aero.Vinf/sqrt(1-AC.Aero.Mach^2).*dy'.*y';
else
    D2 = nan;
    U2 = nan;
    Gamma2 = nan;
    Load2 = nan;
    dA2_dX = nan;
    dRHS2_dX = nan;
    dLtot2_dX = nan;
    dM2_dX = nan;
    dM_dGamma = nan;
    dK2_dx = nan; 
    dF2_dx = nan;
end

%% 2D airfoil drag
if Drag==1
    Zcd = [zeros(1,length(Gamma))'; dCD_dU(:,AC.FEM.ActiveDOFs)'; dCD_dAlpha'; dCD_dAlpha_i'];
    dCL_dX = dCl2d_dX - dClvlm_dX;
    
    % CDi from vlm 
    [dA, dV, ~, dB] = AIC_RHS_Aileron(gradientinit(X),AC,Alpha,U,Xea,DV,DP);
    A = dA.x;
    dA_dX = dA.dx;
    dV_dX = dV.dx;
    dB_dX = dB.dx;
    
   
%     dWind_dX = zeros(size(dV_dX));
    dWind_dX = zeros(nS-1,length(X));
    dGamma_dX= zeros(size(dV_dX));
    
    for i=1:length(X)
        dGamma_dX(:,i) = A\(dV_dX(:,i)-dA_dX(:,:,i)*Gamma);
        dG_dX = [dGamma_dX(end-nS+2:end-1,i)-dGamma_dX(end-nS+3:end,i); dGamma_dX(end,i)]/sqrt(1-AC.Aero.Mach^2);
%         dWind_dX(:,i) = dB_dX(:,:,i)*Gamma + B*dGamma_dX(:,i);
        dWind_dX(:,i) = dB_dX(:,:,i)*G + B*dG_dX;
    end
    
    
    [~, ~, ~, ~, ~, ~, ~, ~, dCDivlm] = Lift_VLM(gradientinit(X),AC,Gamma,Wind,DV,DP);
    pCDivlm_pX = dCDivlm.dx;
    
    dCDivlm_dx = pCDivlm_pX + dCDivlm_dGamma*dGamma_dX + dCDivlm_dWind*dWind_dX;
    pCD_pX = pCD_pX + dCDivlm_dx;
    
else
    Zcd = NaN;
    dCL_dX = nan;
    pCD_pX = nan;
end


%% Matrices

[K, F, M, dK_dx, dF_dx, dM_dx] = Matrices(AC,Load,X,cNode,Xea,Yea,dcNode,dXea,dYea,n,Vf,dVf_dX,Wfload,dWfload_dx,DV,DP);

if Aileron_Eff ==1
    [~, ~, ~, dK2_dx, dF2_dx] = Matrices(AC,Load2,X,cNode,Xea,Yea,dcNode,dXea,dYea,n,Vf,dVf_dX,Wfload,dWfload_dx,DV,DP);
end

%% Stress calculation

[ Failure, ~, Z, dg_dx] =  Stress(AC,U,X,cNode,Xea,Yea,dcNode,dXea,dYea,DV,DP);


%% sensitivity

[Dg_Dx, dU_dX, dMa_dX, dCD_dX] = Adjoint_Sensitivity(AC,X,D,Z,dg_dx,dK_dx,dF_dx,dA_dX,dRHS_dX,dLtot_dX,U,Gamma,dM1_dX, dM2_dX, dM_dGamma,D2,dK2_dx,dF2_dx,dA2_dX,dRHS2_dX,dLtot2_dX,U2,Gamma2, dCL_dX, pCD_pX, Zcd, Drag, Aileron_Eff,dWdes_dx); %dA_da,dAIC_da,d2AIC_dxda, d2V_dxda,si_m);


%% weight

dWw= Weight(AC,gradientinit(X),DV,DP);
dWw_dX = dWw.dx;

dWL = WingLoading(gradientinit(X),AC,DV,DP);
dWL_dX = dWL.dx;

%% Fuel tank volume

dVf = fuel_volume(gradientinit(X),AC,DV,DP);
dVf_dX = dVf.dx;

else
  Dg_Dx = nan;
  dU_dX = nan;
  dMa_dX = nan;
  dCD_dX = nan;  
  dWw_dX = nan;  
  dWL_dX = nan;
  dVf_dX = nan;  
end

