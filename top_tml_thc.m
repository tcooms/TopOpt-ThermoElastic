%% A 167 line code for thermoelastic topology optimization (transient heat conduction) %%
function top_tml_thc(L,h,z,a,Vf,tf,tn,pE,pk,pb,pc)
addpath('C:\Users\...\MMA scripts\'); tic;
%% PARAMETERS & MATERIAL PROPERTIES
[nelx,nely,nele,rmin,dt] = deal(round(L/a),round(h/a),round(L/a)*round(h/a),0.04*L/a,tf/tn);
[E0,E1,nu,alpha,rho0,rho1] = deal(30e-6,30e3,0.3,12e-6,2e-11,2.4e-9);
[k0,k1,cv0,cv1] = deal(0.03,1,rho0*700e6,rho1*900e6);
Fth0 = E1*alpha*z*a/2/(1-nu);
%% PREPARE THERMAL FEA
TA1 = [ 8 -2; -2  8];
TA2 = [-4 -2; -2 -4];
KEth = z/12*[TA1 TA2; TA2 TA1];
TA3 = [4 2; 2 4];
TA4 = [1 2; 2 1];
CE = z/36*a*a*[TA3 TA4; TA4 TA3];
Tnodenrs = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
TedofVec = reshape(Tnodenrs(1:end-1,1:end-1)+1,nele,1);
TedofMat = repmat(TedofVec,1,4)+repmat([0 nely+[1 0] -1],nele,1);
TiK      = reshape(kron(TedofMat,ones(4,1))',16*nele,1);
TjK      = reshape(kron(TedofMat,ones(1,4))',16*nele,1);
Tmaxdof  = (nely+1)*(nelx+1);
% DEFINE LOADS AND SUPPORTS (HEATED BOTTOM EDGE)
Q = zeros(Tmaxdof,1);
T0  = 20; T1 = 800; % [°C]
TT0 = T0*ones(Tmaxdof,1);
T = ones(Tmaxdof,1)*T0;
TBC1 = 1:nely+1:Tmaxdof;
TBC2 = nely+1:nely+1:Tmaxdof-2*(nely+1);
T(TBC2) = T1;
Tc = [TBC1 TBC2];
Tf  = setdiff(1:Tmaxdof,Tc);
Tini = T;
%% PREPARE MECHANICAL FEA
MA11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
MA12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
MB11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
MB12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE  = z/(1-nu^2)/24*([MA11 MA12;MA12' MA11]+nu*[MB11 MB12;MB12' MB11]);
MedofVec = reshape(2*Tnodenrs(1:end-1,1:end-1)+1,nele,1);
MedofMat = repmat(MedofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nele,1);
MiK      = reshape(kron(MedofMat,ones(8,1))',64*nele,1);
MjK      = reshape(kron(MedofMat,ones(1,8))',64*nele,1);
Mmaxdof  = 2*(nely+1)*(nelx+1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
Fm = sparse(2,1,-10000/2,Mmaxdof,1);
F = zeros(Mmaxdof,tn+1);
U = zeros(Mmaxdof,tn+1);
Mc = [Mmaxdof 1:2:nely*2+1];
Mf  = setdiff(1:Mmaxdof,Mc);
%% PREPARE THERMO-MECHANICAL LOAD VECTOR AND ADJOINT VARIABLE
DT = zeros(nele,tn+1);
Fth = zeros(Mmaxdof,tn+1);
Itm = [-1 -1 1 -1 1 1 -1 1];
TT = sparse(repmat((1:nele)',1,4),TedofMat,1/4,nele,Tmaxdof);
TTM = sparse(repmat((1:nele)',1,8),MedofMat,repmat(Itm,nele,1),nele,Mmaxdof);
mu = zeros(Tmaxdof,tn+2);
dc1 = zeros(nele,tn+1);
dc2 = zeros(nele,tn+1);
W = ones(tn+1,1);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
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
xnew  = ones(nely,nelx)*Vf;
xold2 = xnew(:); xold1 = xnew(:); xval  = xnew(:); xmax = ones(nele,1); xmin = zeros(nele,1);
[iter,change,beta,eta,xTilde,upp,low] = deal(0,1,1,0.5,xnew,xmax,xmin);
xnew = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));
%% START ITERATION
while change > 0.01 && iter < 1000
    iter = iter + 1; time0 = toc;
    % MATERIAL INTERPOLATION
    [Ex,Edx] = MaterialInterpolation(xnew,E0,E1,pE);
    [kx,kdx] = MaterialInterpolation(xnew,k0,k1,pk);
    [cx,cdx] = MaterialInterpolation(xnew,cv0,cv1,pc);
    [fx,fdx] = MaterialInterpolation(xnew,E0/E1,1,pb);
    % THERMAL FEA
    TsK = reshape(KEth(:)*kx,16*nele,1);
    Kth = sparse(TiK,TjK,TsK); Kth = (Kth+Kth')/2;
    TsC = reshape(CE(:)*cx,16*nele,1);
    C = sparse(TiK,TjK,TsC); C = (C+C')/2;
    T(:,1) = Tini; Told = Tini; Tnew = Tini; Qnew = zeros(Tmaxdof,1);
    dCK = decomposition(C(Tf,Tf)/dt+Kth(Tf,Tf));
    for i = 1:tn
        Qnew(Tf) = C(Tf,Tf)/dt*Told(Tf)-Kth(Tf,Tc)*Told(Tc)+Q(Tf);
        Tnew(Tf) = dCK\Qnew(Tf);
        Told = Tnew; T(:,i+1) = Tnew;
    end
    % THERMO-MECHANICAL LOAD VECTOR
    DT(:,2:end) = TT*(T(:,2:end)-TT0);
    Fth = TTM'*(Fth0*fx'.*DT);
    % MECHANICAL FEA
    F = Fm + Fth;
    MsK = reshape(KE(:)*Ex,64*nelx*nely,1);
    K = sparse(MiK,MjK,MsK); K = (K+K')/2;
    U(Mf,:) = K(Mf,Mf)\F(Mf,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    v  = mean(xnew(:));
    dv = ones(nele,1)/nele;
    ct = zeros(tn+1,1);
    U0 = U(:,1);
    ce0 = sum((U0(MedofMat)*KE).*U0(MedofMat),2);
    ct(1) = sum(Ex'.*ce0);
    dc1(:,1) = -Edx'.*ce0 + 2*Fth0*(U0(MedofMat)*Itm').*fdx'.*DT(:,1);
    dFdT = (TTM'.*(Fth0*fx))*TT;
    for i=tn+1:-1:2
        Ut_n = U(:,i); Qadj0 = (Ut_n'*dFdT)'; Tt_n = T(:,i); Tt_n1 = T(:,i-1); 
        ce = sum((Ut_n(MedofMat)*KE).*Ut_n(MedofMat),2);
        ct(i) = sum(Ex'.*ce);
        Qadj = (mu(Tf,i+1)'*C(Tf,Tf)/dt)'-2*W(i)*Qadj0(Tf);
        mu(Tf,i) = dCK\Qadj;
        mui = mu(:,i);
        dc1(:,i) = -Edx'.*ce + 2*Fth0*(Ut_n(MedofMat)*Itm').*fdx'.*DT(:,i);
        dc2a = cdx'/dt.*(sum((mui(TedofMat)*CE).*Tt_n(TedofMat),2));
        dc2b = kdx'.*(sum((mui(TedofMat)*KEth).*Tt_n(TedofMat),2));
        dc2c = -cdx'/dt.*(sum((mui(TedofMat)*CE).*Tt_n1(TedofMat),2));
        dc2(:,i) = dc2a + dc2b + dc2c;
    end
    c = sum(W.*ct);
    dc = sum(W'.*dc1+dc2,2);
    if iter == 1; c_0 = c/10; end
    [c,dc] = deal(c/c_0,dc/c_0);
    % FILTERING/PROJECTION OF SENSITIVITIES
    dx = beta*(1-tanh(beta*(xTilde-eta)).^2)/(tanh(beta*eta)+tanh(beta*(1-eta)));
    dc(:) = H*(dc(:).*dx(:)./Hs);
    dv(:) = H*(dv(:).*dx(:)./Hs);
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
    xTilde(:) = (H*xnew(:))./Hs;
    xnew = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));
    if beta < 32 && mod(iter,50) == 0; beta = 2*beta; end
    % PRINT RESULTS
    fprintf([' I: %5i  Obj: %12.2f    Volume: %6.4f    Change: %6.3f    '...
        'IterTime: %7.2f    TotalTime: %9.2f\n'],iter,(c*c_0),v,change,toc-time0,toc);
    % PLOT DESIGN VARIABLES
    imagesc(1-xnew); colormap(gray); caxis([0 1]); axis equal; axis tight; axis off; drawnow;
end
end
function [Mx,Mdx] = MaterialInterpolation(x,M0,M1,pM)
Mx = M0+x(:)'.^pM*(M1-M0);
Mdx = pM*x(:)'.^(pM-1)*(M1-M0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by T. Ooms (May 2023)                       %
%    Department of Structural Engineering and Building Materials           %
%    Ghent University, Technologiepark-Zwijnaarde 60, Zwijnaarde 9052      % 
%    Belgium															   %
%			                                                			   %
% Please sent your comments to: Ticho.Ooms@UGent.be,                       %
%    Wouter.DeCorte@UGent.be, Ruben.VanCoile@UGent.be                      %
%                                                                          %
% The code is intended for educational purposes under GNU General Public   %
% License v3.0 and the theoretical details are discussed in the paper:     %
% Thermoelastic topology optimization of structural components at elevated %
% temperatures considering transient heat conduction                       %
% Online version: https://doi.org/10.1007/s00366-023-01907-7               %
%                                                                          %
% The latest version of this code can be downloaded from the website:      %
% https://github.com/tcooms/TopOpt-ThermoElastic						   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guarantee that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%