function [GtX] = projfullvec(XT, t , grdXT, tf, BufferSize, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% compute G_t(x) = (x - Proj(x-t gradfunx))/t

Ttag = XT.T - t*grdXT.T;
ProjT = buffproj2SE3_strct(Ttag,n,tf, BufferSize );
GtX.L = grdXT.L;
GtX.T = (XT.T - ProjT)/t;

end

