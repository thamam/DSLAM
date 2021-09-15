function [dXTitm1, dXTit] =  dfdXTit(XTitm1,XTit,Ttilde_itm1)


%% for testing 
% XTit        = reshape(Pvecmat(1:16),4,4);
% XTitp1      = reshape(Pvecmat(17:32),4,4);
% Ttilde_it   = reshape(Ttild(1:16),4,4);

%%
Rot = @(M) M(1:3,1:3);
pos  = @(M) M(1:3,4);
poseinv = @(M)[Rot(M).' , -Rot(M).'*pos(M);  zeros(1,3), 1];


T3x3 = [1,0,0,0,0,0,0,0,0; %Transpoxe matrix
        0,0,0,1,0,0,0,0,0;
        0,0,0,0,0,0,1,0,0;
        0,1,0,0,0,0,0,0,0;
        0,0,0,0,1,0,0,0,0;
        0,0,0,0,0,0,0,1,0;
        0,0,1,0,0,0,0,0,0;
        0,0,0,0,0,1,0,0,0;
        0,0,0,0,0,0,0,0,1] ;

%vectorize and express norm(w,fro)^2 as w.'*w
W = eye(4) - XTit*Ttilde_itm1*poseinv(XTitm1);
w = reshape(W(1:3,1:4),[],1); %discard the last row of W

%compute dfdXTitm1
M1 = kron(eye(4),(Rot(XTit)*Rot(Ttilde_itm1)).');
M2 = [ T3x3 , zeros(9,3);
       kron(eye(3),-pos(XTitm1).'), - Rot(XTitm1).'];    
dXTitm1 = -2*M2.'*M1*w; 

%compute dfdXTit
U = kron((Ttilde_itm1*poseinv(XTitm1)), eye(3));
dXTit = -2 * U* w;

% [Oftnrm2, Ogradft] = symslam_ft(XTitm1,XTit,Ttilde_itm1);
end


% function [dXTitm1, dXTit] =  dfdXTit(XTitm1,XTit,Ttilde_itm1)
% 
% 
% %% for testing 
% % XTit        = reshape(Pvecmat(1:16),4,4);
% % XTitp1      = reshape(Pvecmat(17:32),4,4);
% % Ttilde_it   = reshape(Ttild(1:16),4,4);
% 
% %%
% T3x3 = [1,0,0,0,0,0,0,0,0; %Transpoxe matrix
%         0,0,0,1,0,0,0,0,0;
%         0,0,0,0,0,0,1,0,0;
%         0,1,0,0,0,0,0,0,0;
%         0,0,0,0,1,0,0,0,0;
%         0,0,0,0,0,0,0,1,0;
%         0,0,1,0,0,0,0,0,0;
%         0,0,0,0,0,1,0,0,0;
%         0,0,0,0,0,0,0,0,1] ;
% 
% rotmat = @(M) M(1:3,1:3);
% travec  = @(M) M(1:3,4);
% se3vec = @(T) reshape(T(1:3,1:4),[],1);
% poseinv = @(M)[rotmat(M).' , -rotmat(M).'*travec(M);  zeros(1,3), 1];
% 
% %compute dfdXTit
% cvec = se3vec(eye(4) - XTit*Ttilde_itm1*poseinv(XTitm1));
% M1 = kron(eye(4),rotmat(XTit)*rotmat(Ttilde_itm1));
% M2 = [T3x3 , zeros(9,3);
%       kron(eye(3),-travec(XTitm1).'), -rotmat(XTitm1).']; 
% dXTitm1 = -2*M1*M2*cvec; 
% 
% %compute dfdXTitp1
% W1 = kron((Ttilde_itm1*poseinv(XTitm1))', eye(3));
% dXTit = -2 * W1*cvec;
% 
% % [Oftnrm2, Ogradft] = symslam_ft(XTitm1,XTit,Ttilde_itm1);
%  
% 
% end