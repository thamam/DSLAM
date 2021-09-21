function [gX,gXT, gXL, dfdXT , dhdXT] = mexslamgrad(XT,XL, Ttild, Z, uvarrayLeft, uvarrayRight, tf, BufferSize)
%#codegen
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% gx = full gradient of obj
% dfdx gradient of first cost function
% dhdx gradient of second cost function


%% Setup
ts=max(1,(tf-BufferSize));

%prep ind of 1x12 variable of the 4x4 Tit
T_dof = [ ones(3,4); zeros(1,4)];
T12ind = find(T_dof==1);
T16ind = (1:16).';
Ttm1xt12ind = @(t) (t-2)*16 +[T12ind ; T12ind+16];
T3x3 = [1,0,0,0,0,0,0,0,0; %Transpoxe matrix
    0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0,1,0,0;
    0,1,0,0,0,0,0,0,0;
    0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0,1,0;
    0,0,1,0,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0,1] ;

B_T_OPT = [ 1 0  0 0 ;
    0 0 -1 0 ;
    0 1  0 0 ;
    0 0  0 1  ];

OPT_T_B = B_T_OPT;%\eye(4);

P = [eye(2),zeros(2)];
pi = @(s) s/ ([0, 0, 1, 0]*s);

dpids = @(a)  [ 1/a(3), 0,  -a(1)/a(3)^2 , 0;
    0     , 1/a(3) , -a(2)/a(3)^2 , 0;
    0 , 0 , 0, 0;
    0 , 0, -a(4)/a(3)^2, 1/a(3)];

rotmat = @(M) M(1:3,1:3);
travec  = @(M) M(1:3,4);
se3vec = @(T) reshape(T(T12ind),[],1);
poseinv = @(M)[rotmat(M).' , -rotmat(M).'*travec(M);  zeros(1,3), 1];

M     = size(XL,2);
n     = size(XT,2);
gXL   = zeros(size(XL));
dfdXT  = zeros(size(XT));
dhdXT = zeros(size(XT));
dhdXL = zeros(size(XL));

%% compute dfdX (for each f_t
posvecret = @(Tit,t) reshape(Tit((t-1)*16+T16ind),4,4);
% gXT = zeros(tf*16,n);
for i=1:n
    Ti  = XT(:,i);
    Ttilde_i = Ttild(:,i);
    XTitm1 = posvecret(Ti,ts);
    for t=(ts+1):tf
        XTit = posvecret(Ti,t);
        Ttilde_itm1 = posvecret(Ttilde_i,t-1);
        [dXTitm1, dXTit] =  dfdXTit(XTitm1,XTit,Ttilde_itm1);
% %         %%%%% DBG begin
% %         [Oftnrm2, Ogradft] = symslam_ft(XTitm1,XTit,Ttilde_itm1);
% %         dfdTdiff = (Ogradft - [dXTitm1;dXTit]);
% %         dfdTdiffnrm = norm(dfdTdiff);
% %         %%%%% DBG end
        T12updind = Ttm1xt12ind(t);
        %         dfdXT((t-2)*16+(1:32),i) = dfdXT((t-2)*16+(1:32),i) + [dXTitm1;zeros(4,1);dXTit;zeros(4,1)];
        dfdXT(T12updind,i) = dfdXT(T12updind,i) + [dXTitm1;dXTit];
        XTitm1 = XTit; % circ fwd for next iteration
    end
end
%% compute dhdX and dhdL (for each h_t)
for i=1:n
    uviLeft = uvarrayLeft(:,:,i);   % All landmarks measuremetns for i-th agnet
    uviRight = uvarrayRight(:,:,i);   % All landmarks measuremetns for i-th agnet
    zi = Z(:,:,i);     % All lanmdmakrs IDs
    Ti  = XT(:,i);
    for t=ts:tf
        XTit = posvecret(Ti,t);
        XTitInv = poseinv(XTit);
        OTBXTinv = OPT_T_B*XTitInv;
        % since the landmarks numbering starts at 0, we add 1 to all of them
        zid_t = zi(:,t);
        znnznum = nnz(zid_t);
                
        %note that the data here is not indexed by m but just listed in
        % the order the landmarks id appear in zid_t
        uiL_t = uviLeft(:,t);  
        uiR_t = uviRight(:,t);
        dhtdm = zeros(3,znnznum);
        if isempty(zid_t)
            continue
        end
        zid_t_nnz = zid_t(1:znnznum);
        for m =1:znnznum
            mid = zid_t_nnz(m);
            lm = XL(:,mid);
            lmbar = [lm;1];
            %             s = OPT_T_B*(poseinv(XTit)*lmbar);
            s = OTBXTinv*lmbar;
            
            u_itm = [uiL_t(m);uiR_t(m)];
            %compute derivatives using chain rule
            dPPisds = P*dpids(s);
            dhds = 2*(dPPisds).'*(P*pi(s)-u_itm);
            dsdl_m = B_T_OPT(1:4,1:3)*(XTit(1:3,1:3).');
            dhdl_m = dsdl_m.'*dhds ;
            %             dsdT = B_T_OPT(1:4,1:3)*[kron(eye(3),(lm-travec(XTit)).') , -rotmat(XTit).'];
            dsdT = B_T_OPT(1:4,1:3)*[kron(eye(3),(lm-XTit(1:3,4)).') , -XTit(1:3,1:3).'];
            dhdTitm = dsdT.' * dhds;
            %update gradient w.r.t. features
            %             dhdXL(:,mid) = dhdXL(:,mid)+ dhdl_m;
            dhtdm(:,m) = dhdl_m;
            %update gradient w.r.t. poses
            dhdXT(T12ind+(t-1)*16, i) = dhdXT(T12ind+(t-1)*16, i) + dhdTitm;
        end
        dhdXL(:,zid_t_nnz) = dhdXL(:,zid_t_nnz) + dhtdm;
    end
end
gXT = dfdXT + dhdXT;
gXL = dhdXL;
%fix Ti0 to its initial value by zeroing its respective gradients
cancelgradT0 = true;
if cancelgradT0
    gXT(1:16,1:n)=0;
end
gX.T = gXT;
gX.L = gXL;
end


