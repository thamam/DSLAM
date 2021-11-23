%% General Envrionment settings
clear
close all;
clc;
%
addpath(genpath('DslamUtils\'));
addpath(genpath('SymFunc\'));
addpath(genpath('codegen\'));
addpath(genpath('Tukey'));

% SimresultsPath = 'SEKCLIP_NO_NOISE';
SimresultsPath = 'SEKFCLIPPER';
% SimresultsPath = 'SEKFCLIP_lownoise';
COMPUTEXSTAR = true;
vidtag = 'Reg_noise';
featNoiseSigma = 0.5;%how far is L initial guess from actual location

% Control random generator
rng(0,'twister')
BufferSize = 1; %Restriction on howmany many poses to keep updating

RUNMODES = {'singlerun','sequential'};

runmode = RUNMODES{1};
%========================================================%
%% Define usefull handels
rotmat = @(M) M(1:3,1:3);
travec  = @(M) M(1:3,4);
se3vec = @(T) reshape(T(1:3,1:4),[],1);
poseinv = @(M)[rotmat(M).' , -rotmat(M).'*travec(M);  zeros(1,3), 1];

%% Load data
n=3;
[XTstar,tildPveccmat, Zarray, uvarray, XLstar, T, M ] = dataloader(SimresultsPath, n);
Xstar = struct('T',XTstar,'L',XLstar);
% Xstarvec =[XTstar(:);XLstar(:)];
% Lstar = XLstar(:);

%change feat data into array mat array format
SZ = [M, 50]; %The 50 is for max num of features observed by one robot at any given time
[Zout] = convindtoarray(Zarray,SZ);
[ULout,URout] = convmeastoarray(uvarray,SZ);
Uout = struct('left',ULout,'right',URout);

% Remove landmark measurements from T0;
if false
    Zout(:,1,1:n)=0;
end
%% Unroll the measurements incremantally t=1,...,T
% Sim parameters
feathist = [];
goutstat=[];

Lhathist = zeros(3*M,T);
Xposehist = zeros(T*n*16+3*M,T);

Tmat = zeros(T*16,n);
Ttild = [];

maxiter = 500 ;  % max GPA iterations
eta = 1e-6; % stopping criterion
fths = 1e-6;
opt = [maxiter, eta, fths];
%backtracking parameters
alpha = 1;   % choose alpha in (0, 0.5)
beta = 0.5;     % choose beta in (0,1)
gdops = [alpha, beta];
%% Testing data - Computing loss at solution
if true
    % [cost] = dslamcostfun(X0_fsz, tildPveccmat, Zarray, uvarray, n ,T, M);
    %     cost_star = dslamcostfun(Xstarvec, tildPveccmat, Zarray, uvarray, n ,T, M);        
    cost_star = mexcostfunc(Xstar.T, Xstar.L, tildPveccmat, Zout, ULout,URout, n ,T, BufferSize);
    % Compute gradient
    [gX,gXT, gXL, dfdXT , dhdXT] = mexslamgrad(Xstar.T, Xstar.L, tildPveccmat, Zout, ULout,URout,T, BufferSize);
    tf=5;
    opt(1)=500;
    [XstarOut, statOut] = slamPGD(Xstar, tildPveccmat, Zout,Uout, n ,tf, M, BufferSize, opt, gdops);
end
    %% Solve 
if strcmp(runmode,RUNMODES{1}) % singlerunmode
    %% Initialize measurmeents at t==1
    opt(1) = 10000;
    X0.L = XLstar+featNoiseSigma^2*randn(size(XLstar));
    T0   = XTstar(1:16,1:n);   

    X0.T   = computeposetrajectory(T0, tildPveccmat, n);
    %[Xhat, statOut] = slam_mgd(X0, tildPveccmat, Zout,Uout, n ,T, M, T+10, opt, gdops);
    %[Xhat, statOut] = slam_mgd(X0, tildPveccmat, Zout,Uout, n ,T, M, T+10, opt, gdops);
    [Xhat] = slam_mgd(X0, tildPveccmat, Zout,Uout, n ,T, M, T+10, opt, gdops);
    ResArray = {Xhat};
else %run sequential mode
    %% Initialize measurmeents at t==1
    Lhat = zeros(size(XLstar));
    T0 = XTstar(1:16,1:n);
    X0.T = T0;
    X0.L = Lhat;
    
    % Initialize landmarks registry
    newobsfeat = unique(nonzeros(Zout(:,1,1:n)));
    firstobsfeatures = setdiff(newobsfeat,feathist);
    feathist = newobsfeat;
    %Set initial guess for newobserved L as L_m^0 = L_m^* + noise
    Lnewfeat = XLstar(:,firstobsfeatures);
    LnewfeatWithNoise = featNoiseSigma^2*randn(size(Lnewfeat));
    X0.L(:,firstobsfeatures) = Lnewfeat + LnewfeatWithNoise;
    % At time t, we have access to:
    %           1. uv_t and z_t
    %           2. ode from T_itm1 to T_it, i.e.,  Ttild_iTm1
    resltcnt = 0; %counter used to store results
    Xc = X0;
    for tt=2:1:T
        sprintf('New measurement, t= %i \n',tt)
        %% Acquire new measuremetns
        tTildmeas_t = tildPveccmat((1:16)+(tt-2)*16,1:n);
        Ttild = [Ttild ;tTildmeas_t];
        
        %Compute guess for Tit based on Ti{t-1} and Ttild
        for ii=1:n
            Titm1 = reshape(Xc.T((1:16)+(tt-2)*16,ii),4,4);
            Ttildit = reshape(tTildmeas_t(1:end,ii),4,4);
            Xc.T((1:16)+(tt-1)*16,ii) = reshape(Titm1/Ttildit,16,1);
        end
        %% Update landmark registry
        [firstobsfeatures, LnewfeatWithNoise, feathist] = ...
            updatelandmarks(Zout, feathist, XLstar, tt, n,featNoiseSigma);
        Xc.L(:,firstobsfeatures) = LnewfeatWithNoise;
        %% Solve for t=1,2,...
        %[Xhat, statOut(tt,1:6)] = slamPGD(Xc, Ttild, Zout,Uout, n ,tt, M, BufferSize, opt, gdops);
        [Xhat, statOut(tt,1:5)] = slam_mgd(Xc, Ttild, Zout,Uout, n ,tt, M, BufferSize, opt, gdops);
        Xc = Xhat;
        if mod(tt,5)==0
            resltcnt = resltcnt +1;
            ResArray{resltcnt} = Xhat;
        end
    end
    %% Wraps stats
    tblsz = [T 6];
    varTypes = {'double','double','double','double','double','double'};
    statsVarNames = {'itc', 'fXc','grdfXcNrm','GtXcNrm','1/Lk','GXcNorm'};
    statsTbl = array2table(statOut, 'VariableNames',statsVarNames);
end

%%  Use the following code to compare two solutions in SE(3)
% % Compute the perturbed paths
perturbedposetrajec = computeposetrajectory(T0, tildPveccmat, n);

nowstamp = datetime('now','TimeZone','local','Format','d-MMM-y_HH_mm');
makevideo = true;
trajecplotmode.Noise = 1;
trajecplotmode.star = 1;

if makevideo == true
    vidtit = sprintf('simanim_%s_Buf_%i_fs_10_%s',vidtag,BufferSize,nowstamp);
    solanim(ResArray, vidtit, n, BufferSize, perturbedposetrajec,trajecplotmode)  ;
end
saveresults = true;
if saveresults
    save(sprintf('SimResults_%s',nowstamp),'ResArray','Xhat','BufferSize','SimresultsPath')
end
%%Plotting perturbed trajectory

if false
    %%
    h=figure,
    hold all
    plottrajec(perturbedposetrajec,h,'-.');   
    plottrajec(XStar,h,'-')
end
%% Plotting
if false
    robotsColors = {'rs','gs','bs'};
    % Plotting solution estiamtes
    h1 = figure(1);
    clf(h1)
    hold all
    for ti=1:tt
        for ii=1:n
            Tit = reshape(Xhat.T((1:16)+(ti-1)*16,ii),4,4);
            %Plot only translation
            plot3(Tit(1,4), Tit(2,4),Tit(3,4),robotsColors{ii});
        end
    end
    h1.Children.XLim = [-12,12];
    h1.Children.YLim = [-12,12];
end

if false
    robotsColors = {'rs','gs','bs'};
    % Plotting solution estiamtes
    h2 = figure(2);
    clf(h2)
    hold all
    for tt=goutstat(2:end,1)'
        Xposestar_mat = XTstar(1:16*tt,1:3);
        %         Xfeathat_mat = reshape(Lhathist(:,t),3,M);
        for k=goutstat(2,1):tt
            T_i = Xposestar_mat((k-1)*16+(1:16),1:3); % 16x3 matrix [T_1t,T_2t,T_3t]
            for ii=1:n
                Tit = reshape(T_i(:,ii),4,4);
                %Plot only translation
                plot3(Tit(1,4), Tit(2,4),Tit(3,4),robotsColors{ii});
            end
        end
        
    end
    
    h2.Children.XLim = [-12,12];
    h2.Children.YLim = [-12,12];
    title('X^* poses ')
end
