%% A 108-line Matlab code for thermoelastic topology optimization with uniform thermal loads %% 
function top_tml_uniform(L,h,t,z,Vf,rmin,pE,pb)
addpath('C:\Users\tcooms\OneDrive - ugent\C - Onderzoek\2023-01 GCMMA-MMA-code-1.5\GCMMA-MMA-code-1.5'); tic;  
%% PARAMETERS
nelx = round(L/z); nely = round(h/z); nele = nelx*nely;
ini = Vf;
%% MATERIAL PROPERTIES
E0 = 1; Emin = 1e-6;
nu = 0.3;
alpha = 5e-4;
Fth0 = E0*alpha*t*z/2/(1-nu);
%% UNIFORM TEMPERATURE DIFFERENCE
T1 = 10;                
%% PREPARE MECHANICAL FEA
MA11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
MA12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
MB11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
MB12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE  = t/(1-nu^2)/24*([MA11 MA12;MA12' MA11]+nu*[MB11 MB12;MB12' MB11]);
Mnodenrs = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
MedofVec = reshape(2*Mnodenrs(1:end-1,1:end-1)+1,nele,1);
MedofMat = repmat(MedofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nele,1);
MiK      = reshape(kron(MedofMat,ones(8,1))',64*nele,1);
MjK      = reshape(kron(MedofMat,ones(1,8))',64*nele,1);
Mmaxdof  = 2*(nely+1)*(nelx+1);
% DEFINE LOADS AND SUPPORTS (BI-CLAMPED BEAM)
Fm = sparse(2*(nely+1)*ceil((nelx+1)/2),1,-1,Mmaxdof,1);
U = zeros(Mmaxdof,1);
Mcdofs = [1:2*(nely+1) Mmaxdof-2*(nely+1)+1:Mmaxdof];
Mfdofs  = setdiff(1:Mmaxdof,Mcdofs);
%% PREPARE THERMO-MECHANICAL LOAD VECTOR AND ADJOINT VARIABLE
Itm = [-1 -1 1 -1 1 1 -1 1]; 
TTM = sparse(repmat((1:nele)',1,8),MedofMat,repmat(Itm,nele,1),nele,Mmaxdof);
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
    fx = Emin/E0+xnew(:)'.^pb*(1-Emin/E0);
    dfdx = pb*xnew(:)'.^(pb-1)*(1-Emin/E0);
    % THERMO-MECHANICAL LOAD VECTOR
    DT = T1*ones(nele,1);
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
    dc = reshape(dc1+dc2,nely,nelx);
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
    imagesc(1-xnew);
    set(gcf,'position',[150 450 1200 400]);
    colormap(gray); caxis([0 1]); axis equal; axis tight; axis off; drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by T. Ooms (July 2023)                      %
%    Department of Structural Engineering and Building Materials           %
%    Ghent University, Technologiepark-Zwijnaarde 60, Zwijnaarde 9052      % 
%    Belgium															   %
%			                                                			   %
% Please sent your comments to: Ticho.Ooms@UGent.be,                       %
%    Wouter.DeCorte@UGent.be, Ruben.VanCoile@UGent.be                      %
%                                                                          %
% The code is intended for educational purposes under GNU General Public   %
% License v3.0 and the theoretical details are discussed in the paper:     %
% Compliance-based topology optimization of structural components          %
% subjected to thermo-mechanical loading                                   %
% Online version: https://doi.org/10.1007/s00158-023-03563-3               %
%                                                                          %
% This code is an adaptation of the original code for topology             %
% optimization of thermoelastic structures considering steady-state heat   %
% conduction called top_tml_shc.m, where the spatially varying temperature % 
% field is now substituted for a uniform temperature difference.           %
%                                                                          %
% The latest version of this code can be downloaded from the website:      %
% https://github.com/tcooms/TopOpt-ThermoElastic						   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guarantee that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%