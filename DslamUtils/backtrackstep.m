function [Xn, Lk, GLkXn] =  backtrackstep(Xk,f,gk,Ase3,Gse3, opts)                                        
% f - in \R handle to objective function
% F - R^n-> R^n  F=\nabla f
% A_L(x) = Pc( x - (1/L)* \nabla f(x) )
% G_L(x) = L(x - A_L(x)

G_L = @(L) Gse3(Xk, gk, L);
A_L = @(L) Ase3(Xk, gk, L);
if strcmpi(opts.mode,'nonconvex')
    %% The goal is to estimate L_k
    s = min(100,opts.s) ; % 1 s>0 First guess for L_k
    gamma =opts.gamma ;% 0.7  ;% 0<gamma<1
    eta = opts.eta; % 1.5 ;%eta>1
    fk = f(Xk);
    Lk = s;
    GLkXk = G_L(Lk);
    cnt=0;
    while fk - f(A_L(Lk)) < (gamma/Lk)*norm(se32vec(GLkXk))^2
        Lk = eta*Lk;
        GLkXk = G_L(Lk);
        cnt=cnt+1;
    end
    
else %TBD %convex case
end

% Return the projection with the final L_k
Xn = A_L(Lk);
GLkXn = GLkXk;
end