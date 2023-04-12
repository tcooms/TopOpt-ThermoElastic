%% A 138-line Matlab code for topology optimization with steady thermo-mechanical loads %% 
function top_tml_shc(L,h,t,z,Vf,rmin,pE,pk,pb)
addpath('C:\Users\...\MMA scripts\'); tic;  
%% PARAMETERS
nelx = round(L/z); nely = round(h/z); nele = nelx*nely;
ini   = 0.4;
%% MATERIAL PROPERTIES
E0 = 30e3; Emin = 30e-6;
nu = 0.3;
k0 = 1; kmin = 0.03;
alpha = 12e-6;
Fth0 = E0*alpha*t*z/2/(1-nu);
%% PREPARE THERMAL FEA
TA1 = [ 8 -2; -2  8];
TA2 = [-4 -2; -2 -4];
KEth = t/12*[TA1 TA2; TA2 TA1];
Tnodenrs = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
TedofVec = reshape(Tnodenrs(1:end-1,1:end-1)+1,nele,1);
TedofMat = repmat(TedofVec,1,4)+repmat([0 nely+[1 0] -1],nele,1);
TiK      = reshape(kron(TedofMat,ones(4,1))',16*nele,1);
TjK      = reshape(kron(TedofMat,ones(1,4))',16*nele,1);
Tmaxdof  = (nely+1)*(nelx+1);
% DEFINE LOADS AND SUPPORTS (HEATED BOTTOM EDGE)
T0  = 0; T1 = 800;                        
TT0 = T0*ones(Tmaxdof,1);
Q = sparse([],[],0,Tmaxdof,1);         
T = ones(Tmaxdof,1)*T0;            
TBC1 = 1:nely+1:Tmaxdof;
TBC2 = nely+1:nely+1:Tmaxdof;
T(TBC2) = T1;
Tcdofs = [TBC1 TBC2];
Tfdofs  = setdiff(1:Tmaxdof,Tcdofs);
%% PREPARE MECHANICAL FEA
MA11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
MA12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
MB11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
MB12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE  = t/(1-nu^2)/24*([MA11 MA12;MA12' MA11]+nu*[MB11 MB12;MB12' MB11]);
MedofVec = reshape(2*Tnodenrs(1:end-1,1:end-1)+1,nele,1);
MedofMat = repmat(MedofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nele,1);
MiK      = reshape(kron(MedofMat,ones(8,1))',64*nele,1);
MjK      = reshape(kron(MedofMat,ones(1,8))',64*nele,1);
Mmaxdof  = 2*(nely+1)*(nelx+1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
Fm = sparse(2,1,-10000/2,Mmaxdof,1);
U = zeros(Mmaxdof,1);
Mcdofs = [Mmaxdof 1:2:nely*2+1];
Mfdofs  = setdiff(1:Mmaxdof,Mcdofs);
%% PREPARE THERMO-MECHANICAL LOAD VECTOR AND ADJOINT VARIABLE
Itm = [-1 -1 1 -1 1 1 -1 1]; 
TT = sparse(repmat((1:nele)',1,4),TedofMat,1/4,nele,Tmaxdof);
TTM = sparse(repmat((1:nele)',1,8),MedofMat,repmat(Itm,nele,1),nele,Mmaxdof);
mu_adj = zeros(Tmaxdof,1);
%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
i = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        i = i+1;
        iH(i) = e1;
        jH(i) = e2;
        sH(i) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
xnew  = ones(nely,nelx)*ini;
xold2 = xnew(:); xold1 = xnew(:); xval  = xnew(:);
xmax = ones(nele,1); upp = xmax; xmin = zeros(nele,1); low = xmin;
iter = 0; change = 1;
%% START ITERATION
while change > 0.01 && iter < 1000
    iter = iter + 1; time0 = toc;
    % MATERIAL INTERPOLATION
    Ex = Emin+xnew(:)'.^pE*(E0-Emin); 
    dEdx = pE*xnew(:)'.^(pE-1)*(E0-Emin);
    kx = kmin+xnew(:)'.^pk*(k0-kmin);     
    dkdx = pk*xnew(:)'.^(pk-1)*(k0-kmin);
    fx = Emin/E0+xnew(:)'.^pb*(1-Emin/E0);
    dfdx = pb*xnew(:)'.^(pb-1)*(1-Emin/E0);
    % THERMAL FEA
    TsK = reshape(KEth(:)*kx,16*nele,1);
    Kth = sparse(TiK,TjK,TsK); Kth = (Kth+Kth')/2;
    T(Tfdofs) = Kth(Tfdofs,Tfdofs)\(Q(Tfdofs)-Kth(Tfdofs,Tcdofs)*T(Tcdofs));
    % THERMO-MECHANICAL LOAD VECTOR
    DT = TT*(T-TT0);
    Fth = TTM'*(Fth0*fx'.*DT);
    % MECHANICAL FEA
    F = Fm + Fth; 
    MsK = reshape(KE(:)*Ex,64*nele,1); 
    K = sparse(MiK,MjK,MsK); K = (K+K')/2; 
    U(Mfdofs) = K(Mfdofs,Mfdofs)\F(Mfdofs);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    v  = mean(xnew(:));
    dv  = ones(nele,1)/nele;
    ce = sum((U(MedofMat)*KE).*U(MedofMat),2);
    c  = sum(Ex'.*ce);
    dc1 = -dEdx'.*ce;
    dc2 = 2*Fth0*(U(MedofMat)*Itm').*dfdx'.*DT;
    Q_adj = (-2*U'*(TTM'.*(Fth0*fx))*TT)';
    mu_adj(Tfdofs) = Kth(Tfdofs,Tfdofs)\Q_adj(Tfdofs);
    dc3 = dkdx'.*(sum((mu_adj(TedofMat)*KEth).*T(TedofMat),2));
    dc = reshape(dc1+dc2+dc3,nely,nelx);
    if iter == 1; c_scale = c/10; end
    c = c/c_scale;
    dc = dc/c_scale;
    % FILTERING/MODIFICATION OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs); 
    dv(:) = H*(dv(:)./Hs);
    % MMA UPDATE OF DESIGN VARIABLES
    m = 1; n = nele;
    f0val = c; df0dx = dc(:); 
    fval = v/Vf-1; dfdx = dv(:)'/Vf;
    a0 = 1; a1 = zeros(m,1); c1 = 1000*ones(m,1); d1 = ones(m,1);
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a1,c1,d1);
    xold2 = xold1; xold1 = xval; xval = xmma;
    change = max(abs(xval(:)-xold1(:)));
    xnew = reshape(xmma,nely,nelx);
    xnew(:) = (H*xnew(:))./Hs;
    % PRINT RESULTS
    fprintf([' I: %5i  Obj: %12.2f    Volume: %6.4f    Change: %6.3f    '...
        'IterTime: %7.2f    TotalTime: %9.2f\n'],iter,(c*c_scale),v,change,toc-time0,toc);
    % PLOT DESIGN VARIABLES
    figure(1);
    set(gcf,'position',[150 450 1200 400]);
    imagesc(1-xnew);
    colormap(gray); caxis([0 1]); axis equal; axis tight; axis off; drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by T. Ooms (January 2023)                   %
%    Department of Structural Engineering and Building Materials           %
%    Ghent University, Technologiepark-Zwijnaarde 60, Zwijnaarde 9052      % 
%    Belgium															   %
%			                                                			   %
% Please sent your comments to: Ticho.Ooms@UGent.be,                       %
%    Wouter.DeCorte@UGent.be, Ruben.VanCoile@UGent.be                      %
%                                                                          %
% The code is intended for educational purposes under <license details>    %
% and the theoretical details are discussed in the paper:        	       %
% < To be added >                                                          %
%                                                                          %
% The latest version of this code can be downloaded from the website:      %
% < To be added >														   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guarantee that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

