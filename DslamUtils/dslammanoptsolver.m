function [outputArg1,outputArg2] = dslammanoptsolver(inputArg1,inputArg2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Create the problem structure.
manifold = spherefactory(n);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) -x'*(A*x);
problem.egrad = @(x) -2*A*x;      % notice the 'e' in 'egrad' for Euclidean
 
% Numerically check gradient consistency (optional).
checkgradient(problem);
 
% Solve.
[x, xcost, info, options] = trustregions(problem);
 
end

