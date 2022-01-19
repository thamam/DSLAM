% X.T – Tx3 array where X.T((t-1)*(1:16), i) is vec of the i-th robot pose at time t. 
%
% X.L – 3xM array where X.L(:, m) is the (x,y,z)^T coordinates of the m-th landmark. 


function new_pose = mgd_iteration(X, g_eval, funX)

% Line search for step size.
step_size = 0.001;
step_scale = 0.5;
c = 0.1;

% Update landmarks with regular gradient descent.
new_pose.L = X.L - step_size .* g_eval.L;

% Set descent direction for landmarks.
descent_vector.L = -g_eval.L;


% Iterate through columns of poses to handle multiple robots.
for j = 1:size(X.T, 2)

    for i = 1:(length(X.T)/16)

        pose_indices = 1 + (16*(i-1):(16*i - 1));

        new_pose_matrix = calculate_next_pose(reshape(X.T(pose_indices, j), [4 4]),   ...
                                              reshape(g_eval.T(pose_indices), [4 4]), ...
                                              step_size);
        
        new_pose.T(pose_indices, j) = reshape(new_pose_matrix, [16 1]);
        
        descent_direction_matrix = reshape(X.T(pose_indices, j), [4 4]) \ new_pose_matrix;
        descent_vector.T(pose_indices, j) = reshape(descent_direction_matrix, [16 1]);

    end
end

% This needs to be inner product between step direction and gradient.
inner = dot(se32vec(descent_vector), se32vec(g_eval));

% while funX(X) - funX(new_pose) < -step_size * c * inner
%     
% %while funX(X) - funX(new_pose) < 
%     
%     %funX(X) - funX(new_pose)
%     
%     step_size = step_scale * step_size;
%     
%     % Update landmarks with regular gradient descent.
%     new_pose.L = X.L - step_size .* g_eval.L;
% 
%     % Iterate through columns of poses to handle multiple robots.
%     for j = 1:size(X.T, 2)
% 
%         for i = 1:(length(X.T)/16)
% 
%             pose_indices = 1 + (16*(i-1):(16*i - 1));
% 
%             new_pose.T(pose_indices, j) = reshape(calculate_next_pose(reshape(X.T(pose_indices, j), [4 4]),   ...
%                                                                       reshape(g_eval.T(pose_indices), [4 4]), ...
%                                                                       step_size), ...
%                                                [16 1]);
% 
%         end
%     end
% end
    
end


function new_pose = calculate_next_pose(current_pose, g_eval, step_size)

% Test iterate function value with initial step size.
gradient = gradient_projection(current_pose, -step_size*g_eval);
%new_pose = current_pose * exponential_map(-step_size * gradient);
%new_pose = logm(inv(exponential_map(step_size * gradient)) * current_pose);
new_pose =  current_pose * exponential_map(gradient);

end


% Projects the euclidean gradient X to the tangent space at the specified
% point.
function gradient = gradient_projection(manifold_point, e_gradient)

% Projections can be directly computed.
tangent_gradient = [e_gradient(1:3, 4); 
                    (e_gradient(3, 2) - e_gradient(2, 3)) / 2; 
                    (e_gradient(1, 3) - e_gradient(3, 1)) / 2; 
                    (e_gradient(2, 1) - e_gradient(1, 2)) / 2];

% Needs to be mapped to tangent space at desired point - use adjoint
% action.
R = manifold_point(1:3, 1:3);
t = manifold_point(1:3, 4);
adjoint_matrix = [R skew(t) * R;
                  zeros(3) manifold_point(1:3, 1:3)];

gradient = tangent_gradient;
gradient = inv(adjoint_matrix) * tangent_gradient;

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


% function [G_1, G_2, G_3, G_4, G_5, G_6] = create_modified_basis(manifold_point)
% 
% % Definitions for Lie algebra generators (describes basis for space).
% G_1 = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
% G_2 = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
% G_3 = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
% G_4 = [0 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 0];
% G_5 = [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0];
% G_6 = [0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0];
% 
% % Form stacked matrix.
% stacked_basis = [G_1 G_2 G_3 G_4 G_5 G_6];
% 
% % Form adjoint determined by manifold point.
% if size(manifold_point, 1) == size(manifold_point, 2)
%     manifold_point = reshape(manifold_point, [4 4]);
% end
% 
% % Extract relevant pose fields and form adjoint operator.
% R = manifold_point(1:3, 1:3);
% t = manifold_point(1:3, 4);
% adjoint_matrix = [R skew(t) * R;
%                   zeros(3) manifold_point(1:3, 1:3)];
% 
% G_1 = stacked_basis * kron(adjoint_matrix * [1; zeros(5, 1)], eye(4));
% G_2 = stacked_basis * kron(adjoint_matrix * [zeros(1, 1); 1; zeros(4, 1)], eye(4));
% G_3 = stacked_basis * kron(adjoint_matrix * [zeros(2, 1); 1; zeros(3, 1)], eye(4));
% G_4 = stacked_basis * kron(adjoint_matrix * [zeros(3, 1); 1; zeros(2, 1)], eye(4));
% G_5 = stacked_basis * kron(adjoint_matrix * [zeros(4, 1); 1; zeros(1, 1)], eye(4));
% G_6 = stacked_basis * kron(adjoint_matrix * [zeros(5, 1); 1], eye(4));
%               
% end


function skew_matrix = skew(x)

skew_matrix = [0 -x(3) x(2); 
               x(3) 0 -x(1); 
               -x(2) x(1) 0];

end