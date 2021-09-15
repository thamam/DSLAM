function [Xtk] = minoverPose(X0, Ttild, n, tf)

funX = @(X) buffdslamcostfun_strc(X, Ttild, Z, u, n ,tf, M, BufferSize);
Gt = @(XT,t, grdXT) projfullvec(XT, t , grdXT, t, BufferSize, n);
grdfunX = @(XT) buffdslamObjGrad(XT, Ttild, Z, u, tf, BufferSize);

Xc = X0;
maxiter = 5;
fXc = 
[A_Lip] = @(X,Lk) projdlsm_A(X,grdfunX,Lk,Proj ) ;
[G_Lip] = @(X,Lk) projdlsm_G( X,Lk, A_Lip);

while( itrcnt<maxiter G_Lip )
    

    bktrkops.mode = {'nonconvex'};
    [Xn,Lk] =  backtrackstep(Xc,funX,grdfunX, A_Lip, G_Lip, bktrkops);
    
end




end