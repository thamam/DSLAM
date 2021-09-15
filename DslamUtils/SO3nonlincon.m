function [c,ceq] = SO3nonlincon(X,T,n)


%% Parse input variables
T_L_sep = 16*n*T;
Pvecmat = reshape(X(1:T_L_sep),[],3);
%% Computing Poses loss
posvecret = @(ti,t) reshape(ti((t-1)*16+1:t*16),4,4);
cnt=1;
for i=1:n
    Ti  = Pvecmat(:,i);
    for t=1:T        
        Tit = posvecret(Ti,t);
        ceq(cnt,1) = log(det(Tit));
        ceq(cnt+1,1) = norm(Tit*Tit.'-eye(4),'fro');                               
        cnt = cnt+2;
    end   
end
% c = (x(1)-1/3)^2 + (x(2)-1/3)^2 - (1/3)^2;
c=[];