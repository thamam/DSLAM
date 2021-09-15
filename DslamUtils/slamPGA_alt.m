function [Xhat, outstat] = slamPGA_alt(X0, Ttild, Z, u, n ,tf, M, BufferSize, opts, gdops)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if exist('opts','var')
    maxiter = opts(1);
    eta     = opts(2);
    fths    = opts(3);
    %     opt = [maxIter, eta, fths];
else
    eta = 1e-5; % stopping criterion
    fths = 1e-4;
    maxiter = 20 ; % maximum number of GD iterations
end
if exist('gdops','var')
    alpha = gdops(1);
    beta = gdops(2);
else
    alpha = 0.15;   % choose alpha in (0, 0.5)
    beta = 0.5;     % choose beta in (0,1)
end


funX = @(X) buffdslamcostfun_strc(X, Ttild, Z, u, n ,tf, M, BufferSize);



while (itercnt<maxiter ) %alt outer loop
    
    [XTk] = minoverPose
    
    
    [XLk] = minoverLanmrk
    
    % combine and decide if to exit
    
end % alt loop








end

