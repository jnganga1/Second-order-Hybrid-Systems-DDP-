function [params, robot_params,Dyn_params] = set_system_functions(x, u, mode)
addpath([pwd '/support']);
addpath([pwd '/algorithm']);

x_size = length(x);
u_size = length(u);
q_size = length(x)/2;
y_size = 2; % 2D quadruped

params.x_size = x_size;
params.u_size = u_size;
params.q_size = q_size;
params.y_size = y_size;


% Split the state
q = x(1:q_size,1);
qd = x(q_size+1:end,1);

% Construct robot model and obtain H(q) and C(q,qd)
[robot, robot_params] = build_quadruped_model(mode);
[H, C] = HandC(robot, q, qd);
robot_params.x_size = x_size;
robot_params.u_size = u_size;
robot_params.q_size = q_size;
robot_params.y_size = y_size;

% Free Dynamics (flying-phase dynamics)
S = [zeros(q_size-u_size, u_size);
     eye(u_size)]; % selection matrix
tau_b = S*u - C;    % control term, gravity term, centrifugal term and corriolis term
H_x = MatrixPartialD(H, x);
tau_bx = jacobian(tau_b, x);
tau_bu = jacobian(tau_b, u);

% Contact jacobian and its derivative
ft = 1;
bc = 2;
Pos_ctact = fkin_quadruped_2D(q, robot_params);  % forward kinematics (two contact points)
Pos_ctact(2,:)=[];  % eliminate y (2D quadruped in x-z plane)
ft_J = jacobian(Pos_ctact(:,ft),q);
bc_J = jacobian(Pos_ctact(:,bc),q);

ft_Jd = sym(zeros(size(ft_J)));             % jacobian derivates for two contact points
bc_Jd = sym(zeros(size(bc_J)));
for i = 1:size(ft_J,1)
    for j = 1:size(ft_J,2)
        ft_Jd(i,j) = jacobian(ft_J(i,j),q)*qd;
        bc_Jd(i,j) = jacobian(bc_J(i,j),q)*qd;
    end
end

% front stance dynamics
ft_b = [tau_b;ft_Jd*qd];
ft_b_x = jacobian(ft_b,x);
ft_b_u = jacobian(ft_b,u);
ft_K = [H, ft_J';
        -ft_J, zeros(size(ft_J,1),size(ft_J,1))]; % KKT Dynamics matrix
ft_K_x = MatrixPartialD(ft_K, x);                 % tensor
% Dyn_Params.FrontStance.ft_b =


% back stance dynamics
bc_b = [tau_b;bc_Jd*qd];
bc_b_x = jacobian(bc_b,x);
bc_b_u = jacobian(bc_b,u);
bc_K = [H, bc_J';
        -bc_J, zeros(size(bc_J,1),size(bc_J,1))]; % KKT Dynamics matrix
bc_K_x = MatrixPartialD(bc_K, x);                 % tensor

qd_x =  jacobian(qd,x);

% Impact Dynamics
a = [H*qd; zeros(size(ft_Jd*qd))]; % x-
a_x = jacobian(a,x);
a_u = jacobian(a,u); % all zeros
q_x = jacobian(q,x); % unify identity matrix and zeros

% matlabFunction(H, tau_b, H_x, tau_bx, tau_bu, qd_x, 'file', 'support/FreeDynamics', 'vars', {x, u});
% matlabFunction(ft_K,ft_b,ft_K_x,ft_b_x, ft_b_u, qd_x, 'file','support/FrontStanceDyn','vars',{x, u});
% matlabFunction(bc_K,bc_b,bc_K_x,bc_b_x, bc_b_u, qd_x, 'file','support/BackStanceDyn','vars',{x, u});
% 
% matlabFunction(ft_K,a,ft_K_x,a_x, a_u, q_x, 'file','support/FrontImpactDyn','vars',{x, u});
% matlabFunction(bc_K,a,bc_K_x,a_x, a_u, q_x, 'file','support/BackImpactDyn','vars',{x, u});
% matlabFunction(ft_J, ft_Jd,'file','support/FrontJacobian','vars',{x});
% matlabFunction(bc_J, bc_Jd,'file','support/BackJacobian','vars',{x});
end

