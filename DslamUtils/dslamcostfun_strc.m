function [cost] = dslamcostfun_strc(X, Ttild, Zarray, uvarray, n ,T, M)
%Compute loss function for all robots w.r.t. poses and landmarks estimation
% Input:
%   X
% tildPveccmat - ineratial measurements array
% Zarray - landmarks features ids
% uvarray - landmarks features measurements
% n - # of robots
% T - # of poses
% M - # of landmarks
%% Opt variables:
%   * posesarray
%   * Larray
%
% Measurements:
%   * tTarray
%   * Zarray
%% intialize variables
cost = 0;
rotmat = @(M) M(1:3,1:3);
travec  = @(M) M(1:3,4);
se3vec = @(T) reshape(T(1:3,1:4),[],1);
poseinv = @(M)[rotmat(M).' , -rotmat(M).'*travec(M);  zeros(1,3), 1];

%% For debug purpose
% debug = true;
% if debug
%     posehat = posesarray;
% end% of debug
%% Parse input variables
XT = X.T;
XL = X.L;
%% Computing Poses loss
Ftotatlcost=0;
posvecret = @(ti,t) reshape(ti((t-1)*16+1:t*16),4,4);
for i=1:n
    Ti  = XT(:,i);
    Ttilde_i = Ttild(:,i);
    XTit = posvecret(Ti,1);    
    Ki = numel(Ti)/16;
    for t=1:(Ki-1)
        % squeeze(posesarray{1}(t,:,:))  return Tit :4x4 pose at time t
        XTitp1 = posvecret(Ti,t+1);
        % squeeze(tTarray{1}(t,:,:))return tildeTit :4x4 IMU from Tit to Titp1
        %         tildTit = squeeze(tTarray{i}(t,:,:));
        Ttilde_it = posvecret(Ttilde_i,t);
        %         G = Titp1\(tildTit*Tit);
        G = XTitp1*Ttilde_it*poseinv(XTit);
        %         F(t,i) = norm(eye(4) - G,'fro')^2; % Compute the losee
        Ftotatlcost = Ftotatlcost +  norm(eye(4) - G,'fro')^2;
        XTit = XTitp1; % circ fwd for next iteration

    end
end
% Ftotalcost_ = sum(F(:))
%% Computing Landmarks loss
% Define transformations
B_T_OPT = [ 1 0  0 0 ;
    0 0 -1 0 ;
    0 1  0 0 ;
    0 0  0 1 ];

OPT_T_B = B_T_OPT;%\eye(4);

P = [eye(2),zeros(2)];
pi = @(s) s/ ([0, 0, 1, 0]*s);

%compute loss
htotalcost = 0;
hcostarray=zeros(n,T,M);
for i=1:n
    Ki = T;
    uvi = uvarray{i};   % All landmarks measuremetns for i-th agnet
    zi = Zarray{i};     % All lanmdmakrs IDs
    Ti  = XT(:,i);
    for t=1:Ki
        XTit = posvecret(Ti,t);
        % since the landmarks numbering starts at 0, we add 1 to all of them
        zid_t = zi{t};
        %%
        % note that the data here is not indexed by m but just listed in
        % the order the landmarks id appear in zid_t
        zi_t = uvi{t};
        for m =1:numel(zid_t)
            mid = zid_t(m)+1;
            lmbar = [XL(:,mid);1];
            s = OPT_T_B*(XTit\lmbar);
            %             s2 = Tit\(OPT_T_B*lmbar);
            %             s=s2;
            %             s = (Tit\lmbar);
            %             Tit
            %             if abs(s(3))<1e-6
            %                 normalized_l_m = P*s;
            %             else
            normalized_l_m = P*pi(s);
            %             end
            z_itm = zi_t(m,1:2).';
            QL = normalized_l_m - z_itm;
            htotalcost = htotalcost  + norm(QL,2)^2; %h_itm(i,t,zid_t(m)+1);
        end
    end
end
%%
cost = Ftotatlcost + htotalcost;


