function [Pvecmat,tildPveccmat, Zarray, uvarray, L ,T, M] = dataloader(SimresultsPath, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%  IMPORTANT -  Notation and structural conventions     %
%%
% Tvec_K = [ T_10(:); ... ; T_1K;  T_20(:); ... ; T_2K; T_30(:); ... ; T_3K];
% tildTvec_K = [tildTit_10(;) ; ... ; tildTit_1K(;); ... (same as Tvec)
% Lvec_K == Lvec = [l_0(:) ; ... ; l_M(:)]; Lvec size is indepndent of K
%
% X_K = [ Tvec_K ; Lvec_K];
%
%
%
%========================================================%
%========================================================%

%%  Load Distributed Slam simulation data
% T_tilde measurements: advances from T_t to T_{t+1}
inputs = load([SimresultsPath '\inputs.mat']);
% tTarray = {Ttilde1,Ttilde2,Ttilde3};
tildPveccmat=[];
inputsC = struct2cell(inputs);
for ii=1:n
    tildPveci=reshape(shiftdim(inputsC{ii},1),[],1);
    tildPveccmat=[tildPveccmat,tildPveci];
end
% tildPvec_1=reshape(shiftdim(inputs.arr_0,1),[],1);
% tildPvec2=reshape(shiftdim(inputs.arr_1,1),[],1);
% tildPvec3=reshape(shiftdim(inputs.arr_2,1),[],1);
% tildPveccmat = [tildPvec1,tildPvec2,tildPvec3];

% landmarks : Each col is position of landmark in 3D
features_gt = load([SimresultsPath '\features_gt.mat']);
L = features_gt.arr_0;
M = size(L,2);

% IndMask = reshape(1:M*3,3,M); % to be used for landmarks indices


%Poses : each (t,:,:) is the T_t 4x4 pose matrix in SE(3) w.r.t. T_0 = I
posesstruct = load([SimresultsPath '\poses.mat']);
posescarr = struct2cell(posesstruct);
Pvecmat=[];
for ii=1:n
    Pveci=reshape(shiftdim(posescarr{ii},1),[],1);
    Pvecmat=[Pvecmat,Pveci];
end
% Pvec1=reshape(shiftdim(posesstruct.arr_0,1),[],1);
% Pvec2=reshape(shiftdim(posesstruct.arr_1,1),[],1);
% Pvec3=reshape(shiftdim(posesstruct.arr_2,1),[],1);
% Pvecmat = [Pvec1, Pvec2, Pvec3];
% squeeze(posesarray{1}(t,:,:))  return Tit :4x4 pose at time t
%% Visual observations
% id_i_t : array of landmark ids visible for robot i at timestep t . Same
% as the set Vis(L, T_{i,t}). Let M = |Vis(L, T_{i,t})| be the number of
% landmarks visible uv_i_t: 4xM array of landmark observations in the
% stereo camera pair. Same as z_{m,i,t} for all m \in Vis(L, T_{i,t}) . It
% has two z_{m,i,t} one for left camera and one for the right camera in a
% stereo pair. You can just use the left camera (first 2 dims) for
% simplicity.
T = size(posesstruct.arr_0,1);
vis_feat = load([SimresultsPath '\visible_features.mat']);
vis_cell = struct2cell(vis_feat);
U=size(vis_cell,1)/2;  %U= 3*T; U is the index where the last Z_id is stored
uvarray={};
% z_#_id - which landmarks are observed from pose at time t=0,...396
%% Note that landmarks id starts at 0 (not 1 )
for ii=1:n
    Zarray{ii}= vis_cell( ((ii-1)*T+1):ii*T);
    uvarray{ii} = vis_cell(U+(((ii-1)*T+1):ii*T));
end
% z1_id = vis_cell(1:T);
% Z2_id = vis_cell((T+1):2*T);
% z3_id = vis_cell((2*T+1):3*T);
% Zarray = {z1_id, Z2_id, z3_id};

% U = 3*T;
% uv_1 = vis_cell(U+1:U+T);
% uv_2 = vis_cell(U+T+1:U+2*T);
% uv_3 = vis_cell(U+2*T+1:U+3*T);
% uvarray = {uv_1, uv_2, uv_3};


end

