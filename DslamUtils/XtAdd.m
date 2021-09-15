function [XAB] = XtAdd(A,B,tt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
XAB.T = A.T+tt*B.T;
XAB.L = A.L+tt*B.L; 
end

