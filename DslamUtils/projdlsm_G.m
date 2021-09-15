function [G_lip] = projdlsm_G(X,gX,Lk,Pc )
% A_L(x) = Pc( x - (1/L)* \nabla f(x) )  
% G_L(x) = L(x - A_L(x) 

AX = projdlsm_A(X,gX,Lk,Pc);
G_lip.T = Lk*(X.T - AX.T);
G_lip.L = Lk*(X.L - AX.L);

end
