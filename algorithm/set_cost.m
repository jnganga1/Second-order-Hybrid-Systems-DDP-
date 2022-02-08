function set_cost(x, u, y, robot_params,params)

u_size = params.u_size;
x_size = params.x_size;
q_size = params.q_size;
y_size = params.y_size;
q = x(1:q_size,1);

num_mode_gait = params.num_mode_gait;   % get the number of continuous phases
dt = params.dt;
GRF_g = -params.G;   % GRF due to gravity

% state symbolic variables
q = sym('q',[x_size/2 1],'real'); % x, z, pitch, joint angles
qd = sym('qd',[x_size/2 1],'real'); % x, z, pitch, joint angles
x = [q; qd];
y = sym('y',[y_size 1],'real'); % contact force (2x1) (2D quadruped)
u = sym('u',[u_size 1],'real'); % actuations

%% reference 
vd = 0.5;
% bs ref
q_bs_ref = [0,-0.1432,-pi/25,0.35*pi,-0.65*pi,0.35*pi,-0.6*pi]';
qd_bs_ref =[vd;1;zeros(q_size-2,1)]; 
x_bs_ref = [q_bs_ref; qd_bs_ref]; 

% f1 ref
q_f1_ref = [0,-0.1791,pi/35,0.22*pi,-0.65*pi,0.3*pi,-0.7*pi]';
qd_f1_ref = [vd;-1;zeros(q_size-2,1)]; 
x_f1_ref = [q_f1_ref;qd_f1_ref]; 

% fs ref
q_fs_ref = [0,-0.1831,-pi/40,0.35*pi,-0.6*pi,0.33*pi,-0.75*pi]';
qd_fs_ref = [vd;1;zeros(q_size-2,1)]; 
x_fs_ref = [q_fs_ref;qd_fs_ref];          

% f2 ref
q_f2_ref = [0,-0.1785,-pi/25,0.35*pi,-0.7*pi,0.25*pi,-0.65*pi]';
qd_f2_ref =[vd;-1;zeros(q_size-2,1)]; 
x_f2_ref = [q_f2_ref; qd_f2_ref];    

%% set up weighting matrices
% Q weighting on x
% R weighting on u
% S weighting on y

%%%%%%%%%% bs %%%%%%%%
% runnig cost
R_bs = 0.01*diag([5,.5,0.1,0.1]); 
Q_bs = 0.01*diag([0,5,5,2,2,5,5,...
                  8,5,1,20,20,0,0]);
S_bs = 0.1*eye(y_size);      

% terminal cost
Qf_bs = 20*diag([0,.5,5,5,5,5,5,...
                  10,8,0.1,8,8,0,0]);   
              
%%%%%%%%%%% f1 %%%%%%%%
% Phase Running cost
R_f1 = 0.01*eye(u_size);         
Q_f1 = 0.01*diag([0,.5,5,5,5,5,5,...
                  8,2,1,10,10,10,10]); 

% Phase Final cost
Qf_f1 = 5*diag([0,5,8,5,5,5,5,...
                 10,1,0.1,8,8,8,8]); 
                 
%%%%%%%%%% fs %%%%%%%%
% Phase Running cost
R_fs = 0.01*diag([0.1,0.1,5,5]); 
Q_fs= 0.01*diag([0,.5,5,5,5,2,2,...
                 8,5,1,0,0,20,20]);
S_fs = 0.1*eye(y_size);

% Phase Final Cost
Qf_fs = 10*diag([0,.5,5,5,5,5,5,...
                  10,8,0.1,0,0,8,8]);   


%%%%%%%%%%% f2 %%%%%%%%
% Phase Running cost
R_f2 = 0.01*eye(u_size);         
Q_f2 = 0.01*diag([0,.5,5,5,5,5,5,...
                  8,2,1,10,10,10,10]);

% Phase Final cost
Qf_f2 = 5*diag([0,.5,8,5,5,5,5,...
                 10,1,0.1,8,8,8,8]); 



       
%% set up cost
% running cost
l_continuous = ...
    [(x - x_bs_ref)'*Q_bs*(x - x_bs_ref)+u'*R_bs*u+(y - GRF_g)'*S_bs*(y - GRF_g);   % bs
     (x - x_f1_ref)'*Q_f1*(x - x_f1_ref) + u'*R_f1*u;                               % f1 
     (x - x_fs_ref)'*Q_fs*(x - x_fs_ref)+u'*R_fs*u+(y - GRF_g)'*S_fs*(y - GRF_g);   % fs
     (x - x_f2_ref)'*Q_f2*(x - x_f2_ref) + u'*R_f2*u];                              % f2 

l = l_continuous*dt; 

% Augmented Lagrangian
mu = sym('mu',[1,2],'real'); % mu(1) Lagragian multiplier, mu(2) penalty coefficient

syms h1 h2 real % height function of swing foot
h1 = violation_meas(q,robot_params, 2); % front leg touch down
h2 = violation_meas(q,robot_params, 4); % back leg touch down

% terminal cost
lf_phase = ...
    [1/2*(x - x_bs_ref)'*Qf_bs*(x - x_bs_ref);                                          % bs
     1/2*(x - x_f1_ref)'*Qf_f1*(x - x_f1_ref) + 50*(mu(1)*h1 + (1/2*mu(2))^2*h1^2);     % f1    
     1/2*(x - x_fs_ref)'*Qf_fs*(x - x_fs_ref);                                          % fs
     1/2*(x - x_f2_ref)'*Qf_f2*(x - x_f2_ref) + 50*(mu(1)*h2 + (1/2*mu(2))^2*h2^2);];   % f2

LF_no_mu = ...
    [1/2*(x - x_bs_ref)'*Qf_bs*(x - x_bs_ref);   % bs
     1/2*(x - x_f1_ref)'*Qf_f1*(x - x_f1_ref);   % f1    
     1/2*(x - x_fs_ref)'*Qf_fs*(x - x_fs_ref);   % fs
     1/2*(x - x_f2_ref)'*Qf_f2*(x - x_f2_ref);];   % f2

        
  

% Jacobian and hessian
lx = jacobian(l,x);
lu = jacobian(l,u);
ly = jacobian(l,y);
lfx_phase  = jacobian(lf_phase,x);
LF_no_mu_x = jacobian(LF_no_mu,x);

lxx = sym('lxx',[x_size, x_size, num_mode_gait], 'real');
luu = sym('luu',[u_size, u_size, num_mode_gait], 'real');
lux = sym('lux',[u_size, x_size, num_mode_gait], 'real');
lyy = sym('lyy',[y_size, y_size, num_mode_gait], 'real');
lfxx_phase = sym('lfpxx',[x_size, x_size, num_mode_gait], 'real');
LF_no_mu_xx = sym('lf_no_mu_xx',[x_size, x_size, num_mode_gait], 'real');

for i = 1:num_mode_gait 
    lxx(:,:,i) = hessian(l(i,1),x);
    luu(:,:,i) = hessian(l(i,1),u);
    lux(:,:,i) = jacobian(lu(i,:),x);
    lyy(:,:,i) = hessian(l(i,1),y);
    lfxx_phase(:,:,i) = hessian(lf_phase(i,1),x);
    LF_no_mu_xx(:,:,i) = hessian(LF_no_mu(i,1),x);
end


matlabFunction(l,  'vars',{x,y,u},'file','support/running_cost');
matlabFunction(lf_phase, 'vars',{x, mu},'file','support/Phasefinal_cost');
matlabFunction(LF_no_mu, 'vars',{x},'file','support/final_cost_no_mu');

matlabFunction(lx,lu,luu,lxx,lux,ly,lyy, 'vars',{x,y,u},'file','support/lInfo');
matlabFunction(lfx_phase', lfxx_phase, 'vars',{x, mu},'file','support/lfPhaseInfo');
matlabFunction(LF_no_mu_x', LF_no_mu_xx, 'vars',{x},'file','support/lfInfo_no_mu');

end