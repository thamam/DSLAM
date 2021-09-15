function [x] = se32vec(X)
%UNTITLED4 Summary of this function goes here
% For inverse use vec2se3
%   Detailed explanation goes here
x = [X.T(:);X.L(:)];    
end

