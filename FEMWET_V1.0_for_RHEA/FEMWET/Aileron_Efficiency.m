function [eta_a, dM_da, dM_da0, si_m, dA_da, dA_ddelta, dRHS_ddelta] = Aileron_Efficiency(X,AC,D,Wdes,Alpha,U,Xea,delta_a,DV,DP)


ModeNum  = AC.Wing.Airfoils.ModeNum;
nT       = AC.Structure.nT;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
nC       = AC.Q3D.Grid.DnChord;
% nDOF     = AC.FEM.nDOF;


rho  = AC.Aero.rho;
Mach = AC.Aero.Mach;
Q    = AC.Aero.Vinf;


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


nS      = sum(AC.Q3D.Grid.DnSpan);
nElm    = nS-1;
nChord  = sum(nC)-1;
nPanel  = nS-1;

b      = X(4*nT+2*ModeNum*nAirfoil+3);
b1 = AC.Wing.kink*b;
b2 = b - b1;

eta_ai = AC.Aileron.Geom(1); 
eta_af = AC.Aileron.Geom(2);

yN1 = linspace(0,eta_ai,AC.Q3D.Grid.DnSpan(1));
yN2 = linspace(eta_ai,eta_af,AC.Q3D.Grid.DnSpan(2)+1);
yN3 = linspace(eta_af,1,AC.Q3D.Grid.DnSpan(3)+1);
yPanel = [yN1 yN2(2:end) yN3(2:end)]*(b1+b2);

dy = yPanel(2:end) - yPanel(1:end-1);
y = 0.5*(yPanel(2:end) + yPanel(1:end-1));

%% dM_da for rigid wing

% Newton iteration
Gamma = zeros(nPanel*nChord,1);
Alpha0 = 0;

for iter=1:10
    [A, V, dV_dAlpha] = AIC_RHS_Aileron(X,AC,Alpha0,zeros(size(U)),Xea,DV,DP);                        
    Ltot = Lift_VLM(X,AC,gradientinit(Gamma),0,DV,DP);         
    Lift = Ltot.x;
    dL_dGamma = full(Ltot.dx);   
    
    Dr = [A              -dV_dAlpha;
          dL_dGamma       zeros(1,1)];
       
    RHS = -[A*Gamma-V; Lift-Wdes];
    deltas = Dr\RHS;

    dGamma = deltas(1:nPanel*nChord, 1);
    dalpha = deltas(end);
    
    Gamma = Gamma + dGamma;
    Alpha0 = Alpha0 + dalpha;
    
    Error(iter) = norm(deltas, inf); 
    if Error(iter) < 1E-9
        break;
    elseif Error(iter) > 1E10
        fprintf (1, 'Wing is diverging!');
        stop
    end
end
if Error(iter) > 1e-9
    error('The Newton iteration is not converged, increase the maximum number of iterations');
end
    
dM_dGamma = zeros(size(Gamma));
dM_dGamma(end-nS+2:end) = rho*Q/sqrt(1-Mach^2).*dy'.*y';

Z = [dM_dGamma; 0];

Adj = Dr'\Z;
si =  Adj(1:nElm*nChord,:); 

[dA, dRHS] = AIC_RHS_Aileron(X,AC,Alpha0,zeros(size(U)),Xea,DV,DP,gradientinit(delta_a));             
dA_ddelta = dA.dx;
dRHS_ddelta = dRHS.dx;

dA_dar = dA_ddelta*Gamma - dRHS_ddelta;

dM_da0 = - si'*dA_dar; 


%% dM_da for elastic wing 

[dA, dRHS] = AIC_RHS_Aileron(X,AC,Alpha,U,Xea,DV,DP,gradientinit(delta_a));             
A = dA.x;
RHS = dRHS.x;
dA_ddelta = dA.dx;
dRHS_ddelta = dRHS.dx;
Gamma = A\RHS;

dM_dGamma = zeros(size(Gamma));
dM_dGamma(end-nS+2:end) = rho*Q/sqrt(1-Mach^2).*dy'.*y';

Z = [dM_dGamma; zeros(length(U)-6,1); 0]; % zeros(AC.Q3D.N2D,1)];

Adj = D'\Z;
si_m =  Adj(1:nElm*nChord,:);  % phi = Adj(nElm*nChord+1:nElm*nChord+nDOF-6,:);   lambda = Adj(nElm*nChord+nDOF-6+1,:);  psi = Adj(nElm*nChord+nDOF-6+2:end);

dA_da = dA_ddelta*Gamma - dRHS_ddelta;

dM_da = - si_m'*dA_da; %  - phi'*dS_da - lambda'*dW_da; % - psi'*dCL_da; 


%% 

eta_a = dM_da/dM_da0;




