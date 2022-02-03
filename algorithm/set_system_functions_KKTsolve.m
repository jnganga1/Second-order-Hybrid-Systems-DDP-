function [params, robot_params,Dyn_funcs] = set_system_functions_KKTsolve(x, u, mode)
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

import casadi.*
[H, C] = HandC_casadi(robot, q, qd);


robot_params.x_size = x_size;
robot_params.u_size = u_size;
robot_params.q_size = q_size;
robot_params.y_size = y_size;


%% 
[~,~,Dyn_funcs] = AllTogether_set_system_functions_KKTsolve(x,u,mode);
%%


% Free Dynamics (flying-phase dynamics)
S = [zeros(q_size-u_size, u_size);
     eye(u_size)]; % selection matrix
tau_b = S*u - C;    % control term, gravity term, centrifugal term and corriolis term
H_x = jacobian(H,x);

tau_bx = jacobian(tau_b, x); % size 7 by 14
tau_bu = jacobian(tau_b, u); %size 7 by 4

% Contact jacobian and its derivative
ft = 1;
bc = 2;
Pos_ctact = fkin_quadruped_2D(q, robot_params);  % forward kinematics (two contact points)
% Pos_ctact(2,:)=[];
hh = Pos_ctact;
Pos_ctact = [Pos_ctact(1,:); Pos_ctact(3,:)];% eliminate y (2D quadruped in x-z plane)
ft_J = jacobian(Pos_ctact(:,ft),q);
bc_J = jacobian(Pos_ctact(:,bc),q);


ft_Jd = jtimes(ft_J,q,qd); 
bc_Jd = jtimes(bc_J,q,qd);

% front stance dynamics
ft_b = [tau_b;ft_Jd*qd];

%%
ft_b_x = jacobian(ft_b,x);
ft_b_u = jacobian(ft_b,u);
ft_K = [H, ft_J';
        -ft_J, zeros(size(ft_J,1),size(ft_J,1))]; % KKT Dynamics matrix
ft_K_x = jacobian(ft_K,x);%needs to be reshaped

% back stance dynamics
bc_b = [tau_b;bc_Jd*qd];
bc_b_x = jacobian(bc_b,x);
bc_b_u = jacobian(bc_b,u);
bc_K = [H, bc_J';
        -bc_J, zeros(size(bc_J,1),size(bc_J,1))]; % KKT Dynamics matrix
bc_K_x = jacobian(bc_K,x);

qd_x =  jacobian(qd,x);

% Impact Dynamics
a = [H*qd; zeros(size(ft_Jd*qd))]; % x-
a_x = jacobian(a,x);
a_u = jacobian(a,u); % all zeros
q_x = jacobian(q,x); % unify identity matrix and zeros

 
Dyn_funcs.FreeDyn = Function('FreeDyn',{x,u},{H, tau_b, H_x, tau_bx, tau_bu, qd_x});
Dyn_funcs.FrntStncDyn = Function('FrntStncDyn',{x,u},{ft_K,ft_b,ft_K_x,ft_b_x, ft_b_u, qd_x});
Dyn_funcs.BckStncDyn = Function('BckStncDyn',{x,u},{bc_K,bc_b,bc_K_x,bc_b_x, bc_b_u, qd_x});

Dyn_funcs.FrntImpctDyn = Function('FrntImpctDyn',{x,u},{ft_K,a,ft_K_x,a_x, a_u, q_x}); 
Dyn_funcs.BckImpctDyn =  Function('BckImpctDyn',{x,u},{bc_K,a,bc_K_x,a_x, a_u, q_x});
Dyn_funcs.FrntJac = Function('FrntJac',{x},{ft_J, ft_Jd}); 
Dyn_funcs.BckJac = Function('BckJac',{x},{bc_J, bc_Jd});
end
