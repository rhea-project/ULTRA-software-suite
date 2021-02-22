function [K, F, M, dK, dF, dM, dF_dL] = Matrices(AC,Load,Ts,cNode,Xea,Yea,dcNode,dXea,dYea,n,Vf,dVf,Wf,dWf,DV,DP)

global Adjoint


nT       = AC.Structure.nT;
yT       = AC.Structure.yT;
nAirfoil = AC.Wing.Airfoils.nAirfoil;
yAirfoil = AC.Wing.Airfoils.yAirfoil;
nNode    = AC.FEM.nNodea;
nDOF     = AC.FEM.nDOF ;
nElm     = AC.FEM.nElm;
Conc     = AC.FEM.Conc;
nC       = AC.Q3D.Grid.nChord;


%%

K = zeros(nDOF);
F = zeros(nDOF,1);
% M = zeros(nDOF);

dK = zeros(nDOF,nDOF,length(Ts));
dF = zeros(nDOF,length(Ts));
dF_dL = zeros(nDOF,nElm*nC);
% dM = zeros(nDOF,nDOF,length(Ts));

M = nan;
dM = nan;

Gvar = [Ts; Xea; Yea; cNode; Vf; Wf];

for i=1:nElm
  
   if Adjoint ==1
       [Dk,Df] = Element_int(gradientinit(Gvar),AC,Load,nNode,Conc,i,nT,yT,nAirfoil,yAirfoil,nC,n,length(Ts),length(cNode),DV,DP);
       dkv = Dk.dx;
       dfv = Df.dx;
       ke = Dk.x;
       fe = Df.x;
       
       dkt = dkv(:,:,1:length(Ts));
       dkx = dkv(:,:,length(Ts)+1);
       dky = dkv(:,:,length(Ts)+2);
       dkc = dkv(:,:,length(Ts)+3:length(Ts)+2+length(cNode));
       dkf = dkv(:,:,length(Ts)+3+length(cNode));
       dkw = dkv(:,:,length(Ts)+4+length(cNode));
       
       dft = dfv(:,1:length(Ts));
       dfx = dfv(:,length(Ts)+1);
       dfy = dfv(:,length(Ts)+2);
       dfc = dfv(:,length(Ts)+3:length(Ts)+2+length(cNode));
       dff = dfv(:,length(Ts)+3+length(cNode));
       dfw = dfv(:,length(Ts)+4+length(cNode));
  

%        [~,fl] = Element_int(AC,gradientinit(Load),nNode,Conc,i,Ts,nT,yT,nAirfoil,yAirfoil,cNode,Xea,Yea,n);
%        dfl = fl.dx;
        dfl = 0;
     
        
       for t = 1:length(Ts)
           dk(:,:,t) = zeros(12,12);
           df(:,t)   = zeros(12,1); 
           for c = 1:length(cNode)
               dk(:,:,t) = dk(:,:,t)+dkc(:,:,c)*dcNode(c,t);
               df(:,t) = df(:,t) + dfc(:,c)*dcNode(c,t);
           end
%            for l=1:length(Load)
%               df(:,t) = df(:,t) + dfl(:,l)*dL_dX(l,t); 
%            end
          dk(:,:,t) =  dk(:,:,t) + dkt(:,:,t) + dkx*dXea(t) + dky*dYea(t) + dkf*dVf(t) + dkw*dWf(t);
          df(:,t) = df(:,t) + dft(:,t) + dfx*dXea(t) + dfy*dYea(t) + dff*dVf(t) + dfw*dWf(t);
       end
             
%         ke = ky.x;
%         fe = fy;      

   else
%        [ke,fe] = Element_int(AC,Load,nNode,Conc,i,Ts,nT,yT,nAirfoil,yAirfoil,cNode,Xea,Yea,nC);
%        dfl = 0;
       [kl,fl] = Element_int(Gvar,AC,gradientinit(Load),nNode,Conc,i,nT,yT,nAirfoil,yAirfoil,nC,n,length(Ts),length(cNode),DV,DP);
       ke = kl;
       fe = fl.x;
       dfl = fl.dx;
   end
    

    % system matrices
    ii=1+(Conc(i,1)-1)*6;
    jj=1+(Conc(i,2)-1)*6;
    vec1=[ii ii+1 ii+2 ii+3 ii+4 ii+5 jj jj+1 jj+2 jj+3 jj+4 jj+5];
    K(vec1,vec1)= K(vec1,vec1)+ke;
    F(vec1,1)   = F(vec1,1)+fe;
%     M(vec1,vec1)= M(vec1,vec1)+me;

    dF_dL(vec1,:)  = dF_dL(vec1,:) + dfl;
    
    if Adjoint ==1
        dK(vec1,vec1,:)= dK(vec1,vec1,:)+dk;
        dF(vec1,:)     = dF(vec1,:)+df;
%         dM(vec1,vec1,:)= dM(vec1,vec1,:)+dm;
        
    end
end

end



