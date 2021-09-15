function [A_lip] = projdlsm_A(X,gX,Lk,Pc )
% A_L(x) = Pc( x - (1/L)* \nabla f(x) )  
% G_L(x) = L(x - A_L(x) 

T=X.T;
L=X.L;
A_lip.T = Pc(T - (1/Lk)*gX.T);
A_lip.L =   (L - (1/Lk)*gX.L);

end