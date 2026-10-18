function Q_c = computeCoriolisForce(nodes, origin_ref, Psi, DOF, m_i, weights, Z)
n = 7;
Q_c = zeros(n, 1);
A_c = computeRotationMatrix(Z);
phi   = Z(4);
theta = Z(5);
a     = Z(7);
euler_rates = Z(11:13);
dphi   = euler_rates(1);  
dtheta = euler_rates(2);
da     = Z(14);
G = [1,         0,         -sin(theta); 
    0,  cos(phi),  sin(phi)*cos(theta); 
    0, -sin(phi),  cos(phi)*cos(theta)];
omega_prime = G * euler_rates;
G_dot = [0, 0, -cos(theta)*dtheta; 
    0, -sin(phi)*dphi, cos(phi)*cos(theta)*dphi - sin(phi)*sin(theta)*dtheta; 
    0, -cos(phi)*dphi, -sin(phi)*cos(theta)*dphi - cos(phi)*sin(theta)*dtheta];
G_dot_theta_dot = G_dot * euler_rates;
tilde = @(v) [    0, -v(3),  v(2); ...
    v(3),     0, -v(1); ...
    -v(2),  v(1),     0];
w_tilde = tilde(omega_prime);
Gdt_tilde = tilde(G_dot_theta_dot); 
numNodes = size(nodes, 2);
g_vec = [0; 0; -9.80665]; 
[~, idx1] = ismember(1:numNodes, DOF);
[~, idx2] = ismember((1:numNodes) + numNodes, DOF);
[~, idx3] = ismember((1:numNodes) + 2*numNodes, DOF);
Psi_mat = [Psi(idx1).'; Psi(idx2).'; Psi(idx3).']; 
U = (nodes - origin_ref(:)) * 1e-3 + Psi_mat * a; 
U_dot = Psi_mat * da;
w = m_i * weights(:)'; 
A_quad = (w_tilde * w_tilde) * U + 2 * w_tilde * U_dot + Gdt_tilde * U;
F_g_local_unweighted = A_c' * g_vec; 
F_eff = (F_g_local_unweighted - A_quad) .* w; 
Q_c(1:3) = -A_c * (A_quad * w') + g_vec * sum(w);
Q_c(4:6) = G.' * sum(cross(U, F_eff, 1), 2);
Q_c(7) = sum(Psi_mat .* F_eff, 'all'); 
end

% function Q_c = computeCoriolisForce(nodes, origin_ref, Psi, DOF, m_i, weights, Z)
%     n = 7;
%     Q_c = zeros(n, 1);
%     A_c = computeRotationMatrix(Z);
%     % 1. Extract angles, modal coordinates, and their rates
%     phi   = Z(4);
%     theta = Z(5);
%     a     = Z(7);
% 
%     euler_rates = Z(11:13);
%     dphi   = euler_rates(1);  % Correctly extracted from Euler rates, not omega_prime
%     dtheta = euler_rates(2);
%     da     = Z(14);
% 
%     % 2. Build the G Matrix for Z-Y-X (Rz*Ry*Rx)
%     G = [1,         0,         -sin(theta); 
%         0,  cos(phi),  sin(phi)*cos(theta); 
%         0, -sin(phi),  cos(phi)*cos(theta)];
% 
%     % Local angular velocity (omega_prime = G * theta_dot)
%     omega_prime = G * euler_rates;
% 
%     % 3. Calculate time derivative of G (G_dot)
%     G_dot = [0, 0, -cos(theta)*dtheta; 
%              0, -sin(phi)*dphi, cos(phi)*cos(theta)*dphi - sin(phi)*sin(theta)*dtheta; 
%              0, -cos(phi)*dphi, -sin(phi)*cos(theta)*dphi - cos(phi)*sin(theta)*dtheta];
% 
%     % Extra mapping acceleration (G_dot * theta_dot)
%     G_dot_theta_dot = G_dot * euler_rates;
% 
%     % Skew-symmetric matrix operator
%     tilde = @(v) [    0, -v(3),  v(2); ...
%         v(3),     0, -v(1); ...
%         -v(2),  v(1),     0];
% 
%     w_tilde = tilde(omega_prime);
%     numNodes = size(nodes, 2);
%     g_vec = [0; 0; -9.80665]; % Define global gravity vector
% 
%         for i = 1:numNodes
%             Psi_i = [Psi(DOF == i); 
%                 Psi(DOF == i + numNodes); 
%                 Psi(DOF == i + 2 * numNodes)];
% 
%             s_i = (nodes(:,i)' - origin_ref)*1e-3;
%             u_prime = s_i + Psi_i' * a; 
%             u_prime_dot = Psi_i' * da;
% 
%             % a_quad = w~ w~ r_f + 2 w~ c_f_dot - r_f~ G_dot theta_dot
%             a_quad_local = w_tilde * w_tilde * u_prime' + ...
%                        2 * w_tilde * u_prime_dot' - ...
%                        tilde(u_prime) * G_dot_theta_dot;
%             % a_quad_local = w_tilde * w_tilde * u_prime' + 2 * w_tilde * u_prime_dot';
% 
% 
%             % Nodal gravity force
%             F_g_global = m_i * weights(i) * g_vec;
%             F_g_local  = A_c' * F_g_global; % Transform to local frame for moments and modal projections
% 
%             % Accumulate Forces
%             Q_c(1:3) = Q_c(1:3) - m_i * weights(i) * A_c * a_quad_local + F_g_global;           
%             Q_c(4:6) = Q_c(4:6) - m_i * weights(i) * G.' * tilde(u_prime) * a_quad_local + G.' * tilde(u_prime) * F_g_local;  
%             Q_c(7)   = Q_c(7)   - m_i * weights(i) * Psi_i' * a_quad_local + Psi_i' * F_g_local;
%             Q_c = Q_c;
%     end
% end