function [Xsol, it_hist, ierr] = usenewtonsolver(Xc, Ttild, Z,ULout,URout, n ,tf, BufferSize)

%%
tol=[1.d-6, 1.d-6];
% params=[40, 1, 0];
% maxit
% maxdim
% parms = [maxit, maxdim]
% mexslamgrad_mex(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);
Proj = @(XT) mexproj2SE3_strct_mex(XT,n,tf, BufferSize );
grdfunX = @(X) mexslamgrad_mex(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);

vec2se3hand = @(x) vec2se3(x,n,tf);
% Input/outputs needs to be vector 
br_G = @(x) se32vec(projdlsm_G(vec2se3hand(x),grdfunX(vec2se3hand(x)),1,Proj));

br_x = @(X) [X.T(:);X.L(:)];

% br_f = @(x) br_x(mexslamgrad_mex(reshape(x(1:(16*n*tf)),...
%     [],n),reshape(x(1+(16*n*tf):end),3,[]),Ttild, ...
%     Z, ULout,URout, tf, BufferSize));

%%
% [xsol, it_hist, ierr] = brsola(br_x(Xc),br_G,  tol);
[xsol, it_hist, ierr] = nsola(br_x(Xc),br_G,  tol);

XT = reshape(xsol(1:(16*n*tf)),[],n);
XL = reshape(xsol(1+(16*n*tf):end),3,[]);
Xsol = struct('T',XT,'L',XL);
end