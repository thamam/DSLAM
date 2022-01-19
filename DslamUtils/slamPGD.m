function [Xhat, TotOutStat] = slamPGD(X0, Ttild, Z, u, n ,tf, M, BufferSize, opts, gdops)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% ====================================================================%
%%  Projected gradient descent :  find Xhat_t
%  ====================================================================%

%% Parse inputs
[ULout,URout] = deal(u.left,u.right);

%% PGD Algortithm parameters
if exist('opts','var') % opt = [maxIter, eta, fths];
    [maxIter, eta, fths] = deal(opts(1), opts(2), opts(3));
else
    [maxIter, eta, fths] = deal(20, 1e-5, 1e-4);
end

% backtracaking line search parameters: alpha-expected decrease as fraction
% of f(Xc), beta-reduction factor

if exist('gdops','var')
    [alpha, beta] = deal(gdops(1),gdops(2));
else % choose alpha in (0, 0.5),  beta in (0,1)
    [alpha, beta] = deal(0.15, 0.5);    
end

%% Define PGD function handles
%%%function [cost] = mexcostfunc(XT,XL,Ttild,Zarray, uvarrayLeft, uvarrayRight,n,tf,BufferSize)
%funX = @(X) mexcostfunc_mex(X.T,X.L, Ttild, Z, ULout,URout,n ,tf, BufferSize);
funX = @(X) mexcostfunc(X.T,X.L, Ttild, Z, ULout,URout,n ,tf, BufferSize);
%%%function [gX,gXT, gXL, dfdXT , dhdXT] = mexslamgrad(XT,XL, Ttild, Z, uvarrayLeft, uvarrayRight, tf, BufferSize)
grdfunX = @(X) mexslamgrad_mex(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);
%grdfunX = @(X) mexslamgrad(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);
Proj = @(XT) mexproj2SE3_strct_mex(XT,n,tf, BufferSize );
%Proj = @(XT) mexproj2SE3_strct(XT,n,tf, BufferSize );
% Gt = @(XT,t, grdXT) projfullvec(XT, t , grdXT, t, BufferSize, n);


% Set A and G ops w.r.t. the SE3 projection
Ase3 = @(X,grdfX,Lk) projdlsm_A(X, grdfX, Lk, Proj ) ;
Gse3 = @(X,grdfX,Lk) projdlsm_G(X, grdfX, Lk,Proj);
se3vec = @(X) [X.T(:);X.L(:)];    
dist2fp = @(X) norm(se32vec(Gse3(X,grdfunX(X),1)));

%Initialize
Xc = X0;
itc = 0 ;
fXc = funX(Xc);
grdXc = grdfunX(Xc); %buffdslamObjGrad(Xc, Ttild, Z, u, t, BufferSize);
grdXcvec = se3vec(grdXc);
grdfXcNrm = norm(grdXcvec) ;

Lk=1;
G1Xc = Gse3(Xc,grdXc,Lk);
GtXcNrm = norm(se3vec(G1Xc));
%Init stats
outstat = [itc fXc grdfXcNrm GtXcNrm, 1/Lk, -1];
  
bktrkops.s = 1;% s>0 First guess for L_k
bktrkops.gamma =0.95  ;% 0<gamma<1
bktrkops.eta = 2;
bktrkops.mode = {'nonconvex'};
%start PGA

% [Xc, it_hist, ierr] = usenewtonsolver(X0, Ttild, Z,ULout,URout, n ,tf, BufferSize);
dtfpXc = dist2fp(Xc);
while (itc < maxIter &&  dtfpXc > eta) %  dist. to fix-point
    itc = itc+1
    %find Xn with backtracking - B1 in Beck's FOM in opt. book sectopn 10.3.3
    [Xn,Lk,GtXn] =  backtrackstep( Xc,funX,grdXc, Ase3, Gse3, bktrkops);
    bktrkops.s = Lk/4; %update s according to emprirical values of Lk
    fXn = funX(Xn);    
        
    % Compute values for next round
    Xc = Xn;
    fXc=fXn
    grdXc = grdfunX(Xc); 
    grdXcvec = se3vec(grdXc);
    grdfXcNrm = norm(grdXcvec) ;
    GtXc=GtXn;          
    GtXcNrm = norm(se3vec(GtXc));
   
    %Update stats
    outstat(itc+1,:) = [itc fXc grdfXcNrm GtXcNrm, 1/Lk, dist2fp(Xc)];
%     figure,
%     semilogy(outstat(:,1), (outstat(:,2)))

 %%   
%     %% DBG begin
%     [gX,gXT, gXL, dfdXT ,dhdXT] = buffdslamObjGrad(Xc, Ttild, Z, u, tf, BufferSize);
%     funX(Xc)
%     [Lgrad, Posegrad , costVal, PoseGrad_ft, posegrad_ht] = symmultiValComp(Xc,Ttild,tf, Z{1},u{1});%Z{1}=Z{i} i=1 when  n=1;
%      DiffPoseGrad = Posegrad-gXT;
%      diffdhdXT = dhdXT - posegrad_ht;
%      diffdhdXTnorm = norm(diffdhdXT);
%      diffdfdXT = [PoseGrad_ft;0] - dfdXT;
%      diffdfdXTnorm = norm(diffdfdXT);
%     posediffnrm = norm(DiffPoseGrad );
%     DiffLgrad = Lgrad(1:3,1:size(gXL,2)) - gXL;
%     Ldiddnrm = norm(DiffLgrad,'fro')
%     normdiffarr = [diffdhdXTnorm,diffdfdXTnorm, Ldiddnrm];
%     if(sum(normdiffarr)>1e-6)
%         normdiffarr
%     end
%    %%% DBG ends    
end
Xhat = Xc;
dtfpXc = dist2fp(Xc);
TotOutStat = [itc fXc grdfXcNrm GtXcNrm, 1/Lk, dtfpXc];
% [Xnw, it_hist, ierr] = usenewtonsolver(Xnwin, Ttild, Z,ULout,URout, n ,tf, BufferSize);
end



