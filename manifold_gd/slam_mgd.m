function [Xhat, TotOutStat] = slam_mgd(X0, Ttild, Z, u, n ,tf, M, BufferSize, opts, gdops)

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
%grdfunX = @(X) mexslamgrad_mex(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);
grdfunX = @(X) mexslamgrad(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);
Proj = @(XT) mexproj2SE3_strct_mex(XT,n,tf, BufferSize );
Proj = @(XT) mexproj2SE3_strct(XT,n,tf, BufferSize );
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
%Init stats
outstat = [itc fXc grdfXcNrm, 1/Lk, -1];

% [Xc, it_hist, ierr] = usenewtonsolver(X0, Ttild, Z,ULout,URout, n ,tf, BufferSize);
dtfpXc = dist2fp(Xc);
while (itc < maxIter &&  dtfpXc > eta) %  dist. to fix-point
    itc = itc+1;
    
    % Manifold gradient descent with fixed step size.
    Xn = mgd_iteration(Xc, grdXc, 0.0001);
    
    fXn = funX(Xn);    
    % Compute values for next round
    Xc = Xn;
    fXc=fXn;
    grdXc = grdfunX(Xc); 
    grdXcvec = se3vec(grdXc);
    grdfXcNrm = norm(grdXcvec) ;      
   
    %Update stats
    outstat(itc+1,:) = [itc fXc grdfXcNrm, 1/Lk, dist2fp(Xc)];
    
end

Xhat = Xc;
dtfpXc = dist2fp(Xc);
TotOutStat = [itc fXc grdfXcNrm, 1/Lk, dtfpXc];

end
