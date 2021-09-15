function [X] = vec2se3(x,n,tf)
%UNTITLED4 Summary of this function goes here
% For inverse use se32vec
%   Detailed explanation goes here
XT = reshape(x(1:(16*n*tf)),[],n);
XL = reshape(x(1+(16*n*tf):end),3,[]);
X = struct('T',XT,'L',XL);

end


