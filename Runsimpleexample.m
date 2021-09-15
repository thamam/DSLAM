clear ;
close all;
clc;
%%
load SimpleExample.mat
% 
% Pvecmat % poses in vec form each pose is [16,1] , use reshape(Pvecmat,4,4,[]) to change into 4x4X... form
% tildPveccmat   % IMU measurments 
% Zarray  % The t-th cell contains the ids of the landmarks the robot see at time t (numbering start from zero)
% uvarray % Match Zid structure, [x,y] coordinate of the corresponding landmark measurment in the camera plane (Z==1) 

n=1; %number of robots
NT = 10; %number of poses
% Set initial guess 
T0 = Pvecmat(1:16,1:n);

%change feat data into array mat array format
SZ = [M, 50]; %The 50 is for max num of features observed by one robot at any given time
[Zout] = convindtoarray(Zarray,SZ);
[ULout,URout] = convmeastoarray(uvarray,SZ);
Uout = struct('left',ULout,'right',URout);





% Compute the perturbed paths
X0.T = computeposetrajectory(T0, tildPveccmat, n);
% Use noisy version of lanmdrks for initial guess
featNoiseSigma = 0.05;
LobsId =  unique(nonzeros(Zout));
UvNoise = featNoiseSigma^2*randn(3,numel(LobsId));
X0.L(1:3,LobsId) = L(:,LobsId) + UvNoise; %add 1 -lanmdrks numbered from 0


%% Solve using PGA
eta = 1e-6; % stopping criterion
fths = 1e-6;
maxiter = 500 ;
opt = [maxiter, eta, fths];
alpha = 1;   % choose alpha in (0, 0.5)
beta = 0.5;     % choose beta in (0,1)
gdops = [alpha, beta];


[Xhat, stats] = slamPGD(X0, tildPveccmat, Zout,Uout, n ,NT, M, 500, opt, gdops);
