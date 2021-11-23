function [Xhat] = manopt_solve(X0, Ttild, Z, u, n ,tf, M, BufferSize, opts, gdops)

[ULout,URout] = deal(u.left,u.right);

%% PGD Algortithm parameters
if exist('opts','var') % opt = [maxIter, eta, fths];
    [maxIter, eta, fths] = deal(opts(1), opts(2), opts(3));
else
    [maxIter, eta, fths] = deal(20, 1e-5, 1e-4);
end


%% Define PGD function handles
%%%function [cost] = mexcostfunc(XT,XL,Ttild,Zarray, uvarrayLeft, uvarrayRight,n,tf,BufferSize)
%funX = @(X) mexcostfunc_mex(X.T,X.L, Ttild, Z, ULout,URout,n ,tf, BufferSize);
funX = @(X) mexcostfunc(X.T,X.L, Ttild, Z, ULout,URout,n ,tf, BufferSize);
%%%function [gX,gXT, gXL, dfdXT , dhdXT] = mexslamgrad(XT,XL, Ttild, Z, uvarrayLeft, uvarrayRight, tf, BufferSize)
grdfunX = @(X) mexslamgrad_mex(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);
%grdfunX = @(X) mexslamgrad(X.T,X.L, Ttild, Z, ULout,URout, tf, BufferSize);
Proj = @(XT) mexproj2SE3_strct_mex(XT,n,tf, BufferSize );
%Proj = @(XT) mexproj2SE3_strct(XT,n,tf, BufferSize );
% Gt = @(XT,t, grdXT) projfullvec(XT, t , grdXT, t, BufferSize, n);

% Create the problem structure.
elements.T = specialeuclideanfactory(3, length(X0.T)/16);
elements.L = euclideanfactory(3, size(X0.L, 2));
manifold = productmanifold(elements);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) funX(convert_manopt_struct(x));
problem.egrad = @(x) deconvert_manopt_struct(grdfunX(convert_manopt_struct(x)));
options.maxiter = 10000;

% Numerically check gradient consistency (optional).
checkgradient(problem);
 
% Solve.
[x, xcost, info, options] = trustregions(problem, deconvert_manopt_struct(X0), options);

% Display some statistics.
figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');

Xhat = convert_manopt_struct(x);

end

function X = convert_manopt_struct(manopt_struct)

% Fill out X.T segment.
X.T = convert_manopt_T(manopt_struct.T);
X.L = convert_manopt_L(manopt_struct.L);

end


function manopt_struct = deconvert_manopt_struct(X)

% Fill out X.T segment.
manopt_struct.T = deconvert_manopt_T(X.T);
manopt_struct.L = deconvert_manopt_L(X.L);

end


function T = convert_manopt_T(manopt_T)

% Find number of poses.
n = size(manopt_T.R, 3);

% Preallocate T.
T = zeros(n * 16, 1);

for i = 1:n
   
    pose_indices = 1 + (16*(i-1):(16*i - 1));
    T(pose_indices) = reshape([manopt_T.R(:, :, i) manopt_T.t(:, i); zeros(1, 3) 1], [16 1]);
    
end

end


function manopt_T = deconvert_manopt_T(T)

% Find number of poses.
n = size(T, 1) / 16;

for i = 1:n
   
    pose_indices = 1 + (16*(i-1):(16*i - 1));
    
    pose_matrix = reshape(T(pose_indices), [4 4]);
    
    manopt_T.R(:, :, i) = pose_matrix(1:3, 1:3);
    manopt_T.t(:, i) = pose_matrix(1:3, 4);
    
end

end


% X.L – 3xM array where X.L(:, m) is the (x,y,z)^T coordinates of the m-th landmark. 
function L = convert_manopt_L(manopt_L)

L = manopt_L;

end

% X.L – 3xM array where X.L(:, m) is the (x,y,z)^T coordinates of the m-th landmark. 
function manopt_L = deconvert_manopt_L(L)

manopt_L = L;

end