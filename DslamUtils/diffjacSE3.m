
function [grad] =diffjacSE3(Xc, f, f0, VecSE3ID, lndmrksId_T, epsnewVal )
% compute a forward difference Jacobian f'(x), return lu factors
%
% uses dirder.m to compute the columns
%
% C. T. Kelley, November 25, 1993
%
% This code comes with no guarantee or warranty of any kind.
%
%
% inputs:
%         x, f = point and function
%		  f0   = f(x), preevaluated
%

% Temp code
x=Xc;

% Prepare reduced difference indices id - 12X1 for each pose and only
% observed features.
actvId = [VecSE3ID(:); 
            lndmrksId_T(:)];

if exist('epsnewVal','var')
    dirder_h = @( zz) dirder(x,zz,f,f0,epsnewVal);
else
    dirder_h = @( zz) dirder(x,zz,f,f0);
end
n=length(x);
grad = zeros(size(Xc));
for j=actvId.'
    zz=zeros(n,1);
    zz(j)=1;
%     grad(j)=dirder(x,zz,f,f0);
    grad(j)= dirder_h(zz);
end

if size(grad,1)==1
    grad = grad.';
end
% [l, u] = lu(jac);
