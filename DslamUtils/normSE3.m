function [val, cnt] = normSE3(X,Y,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cnt=0;
d=4;
d2=16;
val = 0 ;
if (n>1 && min(size(X))==1)
   X= reshape(X,[],n);
   Y= reshape(Y,[],n);   
end

Kt = size(X,1)/d2;

for i=1:n
   Xi = X(:,i);
   Yi = Y(:,i);
   
   for t=1:Kt %If there are several poses to compare
       Xit = reshape(Xi((1:d2)+(t-1)*d2),d,d);
       Yit = reshape(Yi((1:d2)+(t-1)*d2),d,d);                     
       %seperate to pose and translation
       Xit_p  = Xit(1:d-1,1:d-1);
       Yit_p  = Yit(1:d-1,1:d-1);
       
       Xit_t  = Xit(1:d-1,d);
       Yit_t  = Yit(1:d-1,d);
       
        pnorm_i = norm(eye(3) - Xit_p.'*Yit_p,'fro')^2;
        tnorm_i = norm(Xit_t - Yit_t)^2;
        
        val = val +pnorm_i + tnorm_i;
        cnt=cnt+1;
   end      
    
end

end

