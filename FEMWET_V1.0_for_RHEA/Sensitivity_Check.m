clear all
close all
clc

[AC, T, G, P] = Inputs_A320;

DV = 7; 

cd FEMWET
[Wwing, Failure, U, eta_a, CD, WL, Vf, dWw_dX, Dg_Dx, dU_dX, dMa_dX, dCD_dX, dWL_dX, dVf_dX, CDi, CDp, CDf,dM_da,dM_da0]=FEMWET(AC,T,G,P,6,1,1,0,DV,AC.Weight.FW,AC.Weight.MTOW);
%dM_da,dM_da0

dx = 1e-6;
for n= 1:40

    Xf = [T; G; P; AC.Weight.FW; AC.Weight.MTOW];
    Xf(n) = Xf(n) + dx;
    
    [Wwingf, Failuref, Uf, eta_a_f, CD_f, WLf, Vff,~,~,~,~,~,~,~,~,~,~,dM_daf]= ...
        FEMWET(AC,Xf(1:length(T)),Xf(length(T)+1:length(T)+length(G)),Xf(length(T)+length(G)+1:end-2),6,0,1,0,DV,Xf(end-1),Xf(end));  %dM_daf,dM_da0_f
    
    dCD_dX_f(n) = (CD_f-CD)/dx;
    
    dU_dX_f(:,n) = (Uf-U)/dx;
    
    Dg_Dx_f.U.compression(:,:,n) = (Failuref.U.compression - Failure.U.compression)/dx;
    Dg_Dx_f.U.buckling(:,:,n)    = (Failuref.U.buckling - Failure.U.buckling)/dx;
    Dg_Dx_f.L.tension(:,:,n)     = (Failuref.L.tension - Failure.L.tension)/dx;
    Dg_Dx_f.FS.shear(:,:,n)      = (Failuref.FS.shear - Failure.FS.shear)/dx;
    Dg_Dx_f.FS.buckling(:,:,n)   = (Failuref.FS.buckling - Failure.FS.buckling)/dx;
    Dg_Dx_f.RS.shear(:,:,n)      = (Failuref.RS.shear - Failure.RS.shear)/dx;
    Dg_Dx_f.RS.buckling(:,:,n)   = (Failuref.RS.buckling - Failure.RS.buckling)/dx;
   
    dWw_dX_f(n) = (Wwingf-Wwing)/dx;
    dVf_dX_f(n) = (Vff-Vf)/dx;
    dWL_dX_f(n) = (WLf - WL)/dx;
    
    
    dMa_dX_f(n) = (dM_daf-dM_da)/dx;
end

cd ..