function [gX, gXT, gXL] = dslamObjGrad(X, Ttild, Z, u, tf)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% gx = full gradient of obj
% dfdx gradient of first cost function
% dhdx gradient of second cost function


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
    0 0  0 1 ];

OPT_T_B = B_T_OPT;%\eye(4);

P = [eye(2),zeros(2)];
pi = @(s) s/ ([0, 0, 1, 0]*s);

dpids = @(a)  [ 1/a(3), 0,  -a(1)/a(3)^2 , 0;
    0     , 1/a(3) , -a(2)/a(3)^2 , 0;
    0 , 0 , 0, 0;
    0 , 0, -a(4)/a(3)^2, 1/a(3)];

rotmat = @(M) M(1:3,1:3);
travec  = @(M) M(1:3,4);
se3vec = @(T) reshape(T(1:3,1:4),[],1);
poseinv = @(M)[rotmat(M).' , -rotmat(M).'*travec(M);  zeros(1,3), 1];

M     = size(X.L,2);
n     = size(X.T,2);
gXL   = zeros(size(X.L));
dfdXT  = zeros(size(X.T));
dhdXT = zeros(size(X.T));
dhdXL = zeros(size(X.L));
%% Parse input variables
XT = X.T;
XL = X.L;

%% compute dfdX (for each f_t
posvecret = @(ti,t) reshape(ti((t-1)*16+1:t*16),4,4);
% gXT = zeros(tf*16,n);
for i=1:n
    Ti  = XT(:,i);
    Ttilde_i = Ttild(:,i);
    XTitm1 = posvecret(Ti,1);
    for t=2:tf
        XTit = posvecret(Ti,t);
        Ttilde_itm1 = posvecret(Ttilde_i,t-1);
        [dXTitm1, dXTit] =  dfdXTit(XTitm1,XTit,Ttilde_itm1);
        dfdXT((t-2)*16+(1:32),i) = dfdXT((t-2)*16+(1:32),i) + [dXTitm1;zeros(4,1);dXTit;zeros(4,1)];
        XTitm1 = XTit; % circ fwd for next iteration
    end
end

%% compute dhdX and dhdL (for each h_t)
for i=1:n
    uvi = u{i};   % All landmarks measuremetns for i-th agnet
    zi = Z{i};     % All lanmdmakrs IDs
    Ti  = XT(:,i);
    for t=1:tf
        XTit = posvecret(Ti,t);
        % since the landmarks numbering starts at 0, we add 1 to all of them
        zid_t = zi{t};
        %%
        % note that the data here is not indexed by m but just listed in
        % the order the landmarks id appear in zid_t
        uvi_t = uvi{t};
        for m=1:numel(zid_t)
            mid = zid_t(m)+1;
            lm = XL(:,mid);
            lmbar = [lm;1];
            s = OPT_T_B*(XTit\lmbar);
%             s = OPT_T_B*(poseinv(XTit)*lmbar);
            z_itm = uvi_t(m,1:2).';
            
            %compute derivatives using chain rule
            dPPisds = P*dpids(s);
            dhds = 2*(dPPisds).'*(P*pi(s)-z_itm);                       
            dsdl_m = B_T_OPT(1:4,1:3)*(XTit(1:3,1:3).');           
            dhdl_m = dsdl_m.'*dhds ;      
            dsdT = B_T_OPT(1:4,1:3)*[kron(eye(3),(lm-travec(XTit)).') , -rotmat(XTit).'];
            dhdTitm = dsdT.' * dhds;
            
            %update gradient w.r.t. features
            dhdXL(:,mid) = dhdXL(:,mid)+ dhdl_m;
            
            %update gradient w.r.t. poses
            dhdXT((1:12)+(t-1)*16, i) = dhdXT((1:12)+(t-1)*16, i) + dhdTitm;
        end
    end
end


gXT = dfdXT + dhdXT;

gXL = dhdXL;

%fix Ti0 to its initial value by zeroing its respective gradients
gXT(1:16,1:n)=0;
gX = struct('T',gXT,'L',gXL);
end


