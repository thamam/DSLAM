% X.T – Tx3 array where X.T((t-1)*(1:16), i) is vec of the i-th robot pose at time t. 
%
% X.L – 3xM array where X.L(:, m) is the (x,y,z)^T coordinates of the m-th landmark. 


function poses = cc_mgd_iteration(X, g_eval, funX)

%grid_size = 10;
grid_size = 1;
% grid_range = [0.001 0.5];
% grid_step_scale = 0.5;
poses = cell(1, grid_size);

% Remove later - not needed. (debug)
for i = 1:grid_size
    poses{i} = X;
end


% min_function_value = funX(X) + 1;

%while min_function_value >= funX(X)
    
%     step_sizes = linspace(grid_range(1), grid_range(2), grid_size);
    
    for g = 1:grid_size
    
%         step_size = step_sizes(g);
        step_size = 0.005;
        
        % Need to change this if variable number of landmarks per object (could
        % select subset of indices and set poses.L(indices) = ...).
        poses{g}.L = X.L - step_size .* g_eval.L;

        % Iterate through columns of poses to handle multiple robots.
        for j = 1:size(X.T, 2)

            for i = 1:(length(X.T)/16)

                pose_indices = 1 + (16*(i-1):(16*i - 1));

                poses{g}.T(pose_indices, j) = reshape(calculate_next_pose(reshape(X.T(pose_indices, j), [4 4]),   ...
                                                                       reshape(g_eval.T(pose_indices, j), [4 4]), ...
                                                                       step_size), ...
                                                   [16 1]);

            end

        end
        
        % Debug
%         Proj = @(XT) mexproj2SE3_strct(XT,3,5,1);
%         Ase3 = @(X,grdfX,Lk) projdlsm_A(X, grdfX, Lk, Proj);
%         pgd_solution = Ase3(X, g_eval, 1 / step_size);
%         mgd_diff = poses{1}.T - X.T;
%         pgd_diff = pgd_solution.T - X.T;
%         m_diff_vec = mgd_diff(mgd_diff ~= 0);
%         p_diff_vec = pgd_diff(pgd_diff ~= 0);
%         pgd_cost = funX(pgd_solution);
%         mgd_cost = funX(poses{g});
%         solution_difference = poses{g}.T - pgd_solution.T;
%         solution_difference = solution_difference(solution_difference ~= 0);
        
    end

    % Find smallest function value across grid.
%     function_values = cellfun(funX, poses);

%     [min_function_value, pose_index] = min(function_values);
% 
%     grid_range = grid_range * grid_step_scale;
%end

%plot(step_sizes, function_values); hold on;

% final_step_size = step_sizes(pose_index);
poses = poses{1}; 

end


function new_pose = calculate_next_pose(current_pose, g_eval, step_size)

gradient = gradient_projection(current_pose, g_eval);
new_pose = current_pose * exponential_map(-step_size * gradient);

end


% Projects the euclidean gradient X to the tangent space at the specified
% point.
function gradient = gradient_projection(manifold_point, e_gradient)

% Definitions for Lie algebra generators (describes basis for space).
G_1 = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
G_2 = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
G_3 = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
G_4 = [0 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 0];
G_5 = [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0];
G_6 = [0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0];

% Form stacked matrix.
% stacked_basis = [G_1 G_2 G_3 G_4 G_5 G_6];

% Form adjoint determined by manifold point.
if size(manifold_point, 1) == size(manifold_point, 2)
    manifold_point = reshape(manifold_point, [4 4]);
end

% Extract relevant pose fields and form adjoint operator.
% R = manifold_point(1:3, 1:3);
% t = manifold_point(1:3, 4);
% adjoint_matrix = [R skew(t) * R;
%                   zeros(3) manifold_point(1:3, 1:3)];

% Form modified basis.
% G{1} = stacked_basis * kron(adjoint_matrix * [1; zeros(5, 1)], eye(4));
% G{2} = stacked_basis * kron(adjoint_matrix * [zeros(1, 1); 1; zeros(4, 1)], eye(4));
% G{3} = stacked_basis * kron(adjoint_matrix * [zeros(2, 1); 1; zeros(3, 1)], eye(4));
% G{4} = stacked_basis * kron(adjoint_matrix * [zeros(3, 1); 1; zeros(2, 1)], eye(4));
% G{5} = stacked_basis * kron(adjoint_matrix * [zeros(4, 1); 1; zeros(1, 1)], eye(4));
% G{6} = stacked_basis * kron(adjoint_matrix * [zeros(5, 1); 1], eye(4));
% 
% G{1} = manifold_point * G_1 * inv(manifold_point);
% G{2} = manifold_point * G_2 * inv(manifold_point);
% G{3} = manifold_point * G_3 * inv(manifold_point);
% G{4} = manifold_point * G_4 * inv(manifold_point);
% G{5} = manifold_point * G_5 * inv(manifold_point);
% G{6} = manifold_point * G_6 * inv(manifold_point);
% 
% G{1} = G_1;
% G{2} = G_2;
% G{3} = G_3;
% G{4} = G_4;
% G{5} = G_5;
% G{6} = G_6;

G{1} = manifold_point * G_1;
G{2} = manifold_point * G_2;
G{3} = manifold_point * G_3;
G{4} = manifold_point * G_4;
G{5} = manifold_point * G_5;
G{6} = manifold_point * G_6;

gradient = zeros(6, 1);

for i = 1:6
    gradient(i) = trace(e_gradient.' * G{i}) / norm(G{i}, 'fro')^2;
end
% A = [reshape(G{1}, [16 1]) reshape(G{2}, [16 1]) reshape(G{3}, [16 1]) reshape(G{4}, [16 1]) reshape(G{5}, [16 1]) reshape(G{6}, [16 1])];
% gradient = inv(A.' * A) * A.' * e_gradient(:);

end

function manifold_point = exponential_map(tangent_vector)

u = tangent_vector(1:3);
w = tangent_vector(4:6);
w_x = skew(w);
theta = sqrt(w' * w);

if(theta ~= 0)
    A = sin(theta) / theta;
    B = (1 - cos(theta)) / (theta^2);
    C = (1 - A) / (theta^2);
else
    A = 0;
    B = 0;
    C = 0;
end

R = eye(3) + (A * w_x) + (B * (w_x * w_x));
V = eye(3) + B * w_x+ C *(w_x * w_x);
V_p = V * u;

manifold_point = zeros(4);

manifold_point(1:3,1:3) = R;
manifold_point(1:3,4) = V_p;
manifold_point(4,4) = 1;

end


function skew_matrix = skew(x)

skew_matrix = [0 -x(3) x(2); 
               x(3) 0 -x(1); 
               -x(2) x(1) 0];

end