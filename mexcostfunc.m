function [cost] = mexcostfunc(XT,XL,Ttild,Zarray, uvarrayLeft, uvarrayRight,n,tf,BufferSize)
%#codegen

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
% XTinv = zeros(4,4,n,tf);
%% intialize variables
% cost = 0;
rotmat = @(M) M(1:3,1:3);
travec  = @(M) M(1:3,4);
% se3vec = @(T) reshape(T(1:3,1:4),[],1);
poseinv = @(M)[rotmat(M).' , -rotmat(M).'*travec(M);  zeros(1,3), 1];

%% For debug purpose
% debug = true;
% if debug
%     posehat = posesarray;
% end% of debug
ts=max(1,(tf-BufferSize));
%% Computing Poses loss
Ftotatlcost=0;
posvecret = @(ti,t) reshape(ti((t-1)*16+1:t*16),4,4);
for i=1:n
    Ti  = XT(:,i);
    Ttilde_i = Ttild(:,i);
    XTit = posvecret(Ti,ts);    
    for t=ts:1:tf-1
        XTitp1 = posvecret(Ti,t+1);        
        Ttilde_it = posvecret(Ttilde_i,t);    
        XTitInv = [XTit(1:3,1:3).',  -XTit(1:3,1:3).'*XTit(1:3,4);zeros(1,3), 1];
        G = XTitp1*Ttilde_it*XTitInv;
        Ftotatlcost = Ftotatlcost +  norm(eye(4) - G,'fro')^2;
        XTit = XTitp1; % circ fwd for next iteration
%         XTinv(1:4,1:4,i,t) = XTitInv;
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
% pi = @(s) s/ ([0, 0, 1, 0]*s);
%compute loss
htotalcost = 0;
% hcostarray=zeros(n,tf,M);
for i=1:n
    uviLeft = uvarrayLeft(:,:,i);   % All landmarks measuremetns for i-th agnet
    uviRight = uvarrayRight(:,:,i);   % All landmarks measuremetns for i-th agnet
    zi = Zarray(:,:,i);     % All lanmdmakrs IDs
    Ti  = XT(:,i);
    for t=ts:1:tf
        XTit = posvecret(Ti,t);
        XTitInv = poseinv(XTit);
%         XTitInv = XTinv(1:4,1:4,i,t); 
        OTBXTinv = OPT_T_B*XTitInv;
        % since the landmarks numbering starts at 0, we add 1 to all of them
        zid_t = nonzeros(zi(:,t));
        featNum = numel(zid_t);
        %%
        % note that the data here is not indexed by m but just listed in
        % the order the landmarks id appear in zid_t
        uiL_t = uviLeft(1:featNum,t);  
        uiR_t = uviRight(1:featNum,t);          
        for m =1:featNum
            mid = zid_t(m);
            lmbar = [XL(:,mid);1];
            s = OTBXTinv*lmbar;
            normalized_l_m = P*(s/s(3));
            u_itm = [uiL_t(m);uiR_t(m)];
            QL = normalized_l_m - u_itm;
            htotalcost = htotalcost  + norm(QL,2)^2; %h_itm(i,t,zid_t(m)+1);
        end
    end
end
%%
cost = Ftotatlcost + htotalcost;


