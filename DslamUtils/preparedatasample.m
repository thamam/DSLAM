function [] = preparedatasample()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% General Envrionment settings
%
addpath(genpath('DslamUtils\'));
addpath(genpath('SymFunc\'));
addpath(genpath('codegen\'));

SimresultsPath = 'SEKFCLIPPER';


[Pvecmat,tildPveccmat, Zarray, uvarray, L, T, M ] = dataloader(SimresultsPath,1);
%%

Nred = 10 ;% number of desired measurements

Pvecmat = Pvecmat(1:16*Nred);
tildPveccmat   = tildPveccmat(1:16*(Nred-1));
Zarray  = {Zarray{1}(1:Nred)}
uvarray = {uvarray{1}(1:Nred)};

save('SimpleExample.mat')

end

