function [params, robot_params,Dyn_funcs] = AllTogether_set_system_functions_ExtendedRNEA(x, u, mode)
addpath([pwd '/support']);
addpath([pwd '/algorithm']);

x_size = length(x);
u_size = length(u);
q_size = length(x)/2; Nb=q_size;
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

if isa(q,'sym')
    [H, C] = HandC(robot, q, qd);
else 
    import casadi.*
    [H, C] = HandC_casadi(robot, q, qd);
end

robot_params.x_size = x_size;
robot_params.u_size = u_size;
robot_params.q_size = q_size;
robot_params.y_size = y_size;

%%

% Free Dynamics (flying-phase dynamics)
S = [zeros(q_size-u_size, u_size);
     eye(u_size)]; % selection matrix
tau_b = S*u - C;    % control term, gravity term, centrifugal term and corriolis term

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

% back stance dynamics
bc_b = [tau_b;bc_Jd*qd];

% Impact Dynamics
Im_b = [H*qd; zeros(size(ft_Jd*qd))]; % x-

%front Dynamics KKT matrix
ft_K = [H, ft_J';
        -ft_J, zeros(size(ft_J,1),size(ft_J,1))]; % KKT Dynamics matrix
%back Dynamics KKT matrix
bc_K = [H, bc_J';
        -bc_J, zeros(size(bc_J,1),size(bc_J,1))]; % KKT Dynamics matrix

Dyn_funcs.ft_K = Function('ft_K',{q,qd},{ft_K});

%% This is a place holder for the later Hinv Function -> This will soon be replaced with a function call 

aVec = MX.sym('aVec',[9,1]); 
a7Vec = MX.sym('aVec',[7,1]); 

%Explicit formation of Free Hinv Matrix
Hinv = H \ a7Vec; 
fr_Kinv_vec = Function('fr_Hinv_vec',{x,a7Vec},{Hinv});

Hinv_T = a7Vec' / H; 
fr_Kinv_vec_T = Function('fr_Hinv_vec',{x,a7Vec},{Hinv_T});

ft_Kinv = ft_K \ aVec; 
ft_Kinv_vec = Function('ft_Kinv_vec',{x,aVec},{ft_Kinv}); 

ft_Kinv_T =  aVec' / ft_K; 
ft_Kinv_vec_T = Function('ft_Kinv_vec',{x,aVec},{ft_Kinv_T}); 

bc_Kinv = bc_K \ aVec; 
bc_Kinv_vec = Function('bc_Kinv_vec',{x,aVec},{bc_Kinv });

bc_Kinv_T =  aVec' / bc_K; 
bc_Kinv_vec_T = Function('bc_Kinv_vec',{x,aVec},{bc_Kinv_T});

%% First Derivatives

    %Front Stance Dynamics
    % ft_dyn = ft_K\ft_b; 
    ft_dyn = ft_Kinv_vec(x,ft_b); 
    Dyn_funcs.ft_dyn = Function('ft_dyn',{x,u},{ft_dyn});
    ft_dyn_x_w_avec = ft_Kinv_vec(x,jacobian(ft_b,x) - jacobian(ft_K*aVec,x));
    ft_dyn_x_w_avec_fn = Function('ft_dyn_x_w_avec_fn',{x,aVec},{ft_dyn_x_w_avec});
    ft_dyn_x_soln = ft_dyn_x_w_avec_fn(x,ft_dyn);
    Dyn_funcs.ft_dyn_x = Function('ft_dyn_x_soln',{x,u},{ft_dyn_x_soln});

    ft_dyn_u = ft_Kinv_vec(x,jacobian(ft_b,u));
    Dyn_funcs.ft_dyn_u = Function('ft_dyn_u',{x,u},{ft_dyn_u});

    %Back Stance Dynamics
    % bc_dyn = bc_K \bc_b; 
    bc_dyn = bc_Kinv_vec(x,bc_b); 
    Dyn_funcs.bc_dyn = Function('bc_dyn',{x,u},{bc_dyn});
    bc_dyn_x_w_avec = bc_Kinv_vec(x,jacobian(bc_b,x) - jacobian(bc_K*aVec,x));
    bc_dyn_x_w_avec_fn = Function('bc_dyn_x_w_avec_fn',{x,aVec},{bc_dyn_x_w_avec});
    bc_dyn_x_soln = bc_dyn_x_w_avec_fn(x,bc_dyn); 
    Dyn_funcs.bc_dyn_x = Function('bc_dyn_x_soln',{x,u},{bc_dyn_x_soln});

    bc_dyn_u = bc_Kinv_vec(x,jacobian(bc_b,u));
    Dyn_funcs.bc_dyn_u = Function('bc_dyn_u',{x,u},{bc_dyn_u});

    %Free Dynamics
    % fr_dyn = H \ tau_b; 
    fr_dyn = fr_Kinv_vec(x,tau_b); 
    Dyn_funcs.fr_dyn = Function('fr_dyn',{x,u},{fr_dyn});
    fr_dyn_x_w_avec = fr_Kinv_vec(x,jacobian(tau_b,x) - jacobian(H*a7Vec,x));
    fr_dyn_x_w_avec_fn = Function('fr_dyn_x_w_avec_fn',{x,a7Vec},{fr_dyn_x_w_avec});
    fr_dyn_x_soln = fr_dyn_x_w_avec_fn(x,fr_dyn); 
    Dyn_funcs.fr_dyn_x = Function('fr_dyn_x_soln',{x,u},{fr_dyn_x_soln});

    fr_dyn_u = fr_Kinv_vec(x,jacobian(tau_b,u));
    Dyn_funcs.fr_dyn_u = Function('fr_dyn_u',{x,u},{fr_dyn_u});

    %Front Impact Dynamics h
    % ft_Im_Dyn = ft_K\Im_b; 
    ft_Im_Dyn = ft_Kinv_vec(x,Im_b);
    Dyn_funcs.fr_Im_Dyn = Function('fr_Im_Dyn',{x},{ft_Im_Dyn});
    ftIm_dyn_x_w_avec = ft_Kinv_vec(x,jacobian(Im_b,x) - jacobian(ft_K*aVec,x));
    ftIm_dyn_x_w_avec_fn = Function('ftIm_dyn_x_w_avec_fn',{x,aVec},{ftIm_dyn_x_w_avec});
    ftIm_dyn_x_soln = ftIm_dyn_x_w_avec_fn(x,ft_Im_Dyn); 
    Dyn_funcs.ftIm_dyn_x = Function('ft_dyn_x_soln',{x},{ftIm_dyn_x_soln});

    ftIm_dyn_u = ft_Kinv_vec(x,jacobian(Im_b,u));
    Dyn_funcs.ftIm_dyn_u = Function('ft_dyn_u',{x},{ftIm_dyn_u});

    %Back Impact Dynamics 
    % bc_Im_Dyn = bc_K \Im_b; 
    bc_Im_Dyn  = bc_Kinv_vec(x,Im_b);
    Dyn_funcs.bc_Im_Dyn = Function('bc_Im_Dyn',{x},{bc_Im_Dyn});
    bcIm_dyn_x_w_avec = bc_Kinv_vec(x,jacobian(Im_b,x) - jacobian(bc_K*aVec,x));
    bcIm_dyn_x_w_avec_fn = Function('ftIm_dyn_x_w_avec_fn',{x,aVec},{bcIm_dyn_x_w_avec});
    bcIm_dyn_x_soln = bcIm_dyn_x_w_avec_fn(x,bc_Im_Dyn); 
    Dyn_funcs.bcIm_dyn_x = Function('ft_dyn_x_soln',{x},{bcIm_dyn_x_soln});

    bcIm_dyn_u = bc_Kinv_vec(x,jacobian(Im_b,u)); %This is all zeros 
    Dyn_funcs.bcIm_dyn_u = Function('ft_dyn_u',{x},{bcIm_dyn_u});

%%
%% Second Derivatives 
%% The things we need first

%Our Rotations happen in 3d space, adding all zero y position
Pos_ctact_3d = [hh(1,:); 0*hh(2,:) ; hh(3,:)]; 

mu = MX.sym('mu',Nb); % represents qqd part of muHinv
muf = MX.sym('muf',[2 1]); %represents f part of muHinv
        
%The qdd parts of first_qMX
nu_q = MX.sym('nu_q',[Nb,Nb]);
nu_qd = MX.sym('nu_qd',[Nb,Nb]);
nu_tau = MX.sym('nu_tau',[Nb,u_size]);

%The lambda parts of first_qMX
f_q  = MX.sym('f_q',[2,7]);
f_qd  = MX.sym('f_qd',[2,Nb]); 
f_tau = MX.sym('f_tau',[2,u_size]); 

%Zeroing out the y dimension since we are in the 2d plane
f_Q  = [f_q(1,:);MX.zeros(1,7);f_q(2,:)];
f_Qd  = [f_qd(1,:);MX.zeros(1,7);f_qd(2,:)];
f_Tau = [f_tau(1,:);MX.zeros(1,u_size);f_tau(2,:)];

spat_f_q = MX.zeros(6,7);
spat_f_qd = MX.zeros(6,7);
spat_f_tau = MX.zeros(6,u_size);
for idx = 1:Nb 
    %Fpt assumes a (x,y,z) point and an (fx,fy,fz) force. Since our "Force"
    %is a matrix (i.e. partials), we recursively use the columns as forces.
    spat_f_q(:,idx) = Fpt(-f_Q(:,idx),Pos_ctact_3d(:,ft)); %Introduces jacobian to f_q  
    spat_f_qd(:,idx) = Fpt(-f_Qd(:,idx),Pos_ctact_3d(:,ft)); %Introduces jacobian to f_q  
end
for idx = 1:u_size
    spat_f_tau(:,idx) = Fpt(-f_Tau(:,idx),Pos_ctact_3d(:,ft)); %Introduces jacobian to f_q  
end 

%Format Required by RNEA algo ::state 5 is front legs. state 7 is hind legs
ft_q_ext_MX = cell([1 Nb]); ft_q_ext_MX{5} = spat_f_q;
ft_qd_ext_MX = cell([1 Nb]); ft_qd_ext_MX{5} = spat_f_qd;
ft_tau_ext_MX = cell([1 Nb]); ft_tau_ext_MX{5} = spat_f_tau;

bc_q_ext_MX = cell([1 Nb]); bc_q_ext_MX{7} = spat_f_q;
bc_qd_ext_MX = cell([1 Nb]); bc_qd_ext_MX{7} = spat_f_qd;
bc_tau_ext_MX = cell([1 Nb]); bc_tau_ext_MX{7} = spat_f_tau;


%zero out the robot's gravity 
robot_no_grav = robot; robot_no_grav.gravity = 0*robot_no_grav.gravity;
%since f_q is 2D but robot expressed in 3D 
% % [mujacPart_q,~,~,~] = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_qd_ext_MX);
%  mujacPart_q = mu'*ft_J' * f_qd; 
 

  qVal = rand(size(q)); qdVal = rand(size(qd));
  tauVal = rand(size(u)); muVal = rand(size(mu)); 
  xVal = [qVal; qdVal];  
  
  h1 = Dyn_funcs.ft_dyn_x(xVal,tauVal);  h2 = Dyn_funcs.ft_dyn_u(xVal,tauVal);
  all = [h1 h2];   FOPVal = all(1:Nb,:);
  fQVal = all(Nb+1:end,Nb+1:Nb*2);
  ftauVal = h2(Nb+1:end,:);
  
%  Tried = Function('Tried',{q,qd,mu,f_qd},{mujacPart_TRIED});
%  Trued = Function('Trued',{q,qd,mu,f_qd},{mujacPart_q}) 
%  aa = Tried(qVal,qdVal,muVal,fQVal) 
%  bb = Trued(qVal,qdVal,muVal,fQVal) 
%  1==1;

%
%% EXPLICIT METHODS
%% Front Stance Dynamics 
  %Eqn 32 part c
  %wrt q
  mujacPart_q = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_q_ext_MX);
  muMpart_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q);  %This is correct  
  ExtraParts_q  = mujacPart_q + muMpart_q' - muf'*ft_J*nu_q;

  %wrt qd
  mujacPart_qd = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_qd_ext_MX);
  muMpart_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd);  %This is correct  
  ExtraParts_qd  = mujacPart_qd + muMpart_qd' - muf'*ft_J*nu_qd;

  %wrt tau
  mujacPart_tau = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_tau_ext_MX);
  muMpart_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau);  %This is correct  
  ExtraParts_tau  = mujacPart_tau + muMpart_tau' - muf'*ft_J*nu_tau;

  %eqn 32 a 
  mu_rhs_q = jacobian([mu;muf]'*ft_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*ft_b,qd);
  
  %eqn 32 b 
  muHvec = [mu;muf]'*ft_K*aVec; 
  
  %Derivs 
  part31a_Brb = jacobian(mu_rhs_q,q); 
  part31b_Brb = jacobian(jacobian(muHvec,q),q);
  part31c_Brb = jacobian(ExtraParts_q,q); % jtimes([jac*nu_q]',q, muf);% + muf'*jac*nu_q,q); 
  ft_out_qq = part31a_Brb - (part31c_Brb+part31c_Brb') - part31b_Brb;
  
  %This is qqd  - works
  part31a_Brb = jacobian(mu_rhs_q,qd); 
  part31b_Brb = jacobian(jacobian(muHvec,q),qd); % + jacobian(jacobian(muf'*a_df_end,q),qd);
  part31c_Brb = jacobian(ExtraParts_qd,q); %+ muf'*jac*nu_q,q); 
  ft_out_qqd = part31a_Brb - (part31c_Brb)'-part31b_Brb;
  
  ft_out_qdqd = jacobian(mu_rhs_qd,qd); %Eqn 31a  
  ft_out_qtau = -jacobian(ExtraParts_tau,q);%+ muf'*jac*nu_tau,q);
  ft_secondMat  = [ft_out_qq ft_out_qqd; ft_out_qqd' ft_out_qdqd]; 
  
  FOP = [nu_q nu_qd nu_tau];
  f_Q = [f_q f_qd f_tau];
  bb=[mu;muf];
  second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{ft_secondMat,ft_out_qtau}); 
  [ft_OutSecondMat, ft_OutQTau] = second_mat_fn(q,qd,u,ft_dyn,ft_Kinv_vec_T(x,bb),FOP,f_Q); 
  ft_OutSecondMat_fn = Function('OutSecondMat_fn',{q,qd,u,bb,FOP,f_Q},{ft_OutSecondMat,ft_OutQTau}); 

%% Back Stance Dynamics 
  %Eqn 32 part c
  %wrt q
  mujacPart_q = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,bc_q_ext_MX);
  muMpart_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q);  %This is correct  
  ExtraParts_q  = mujacPart_q + muMpart_q' - muf'*bc_J*nu_q;

  %wrt qd
  mujacPart_qd = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,bc_qd_ext_MX);
  muMpart_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd);  %This is correct  
  ExtraParts_qd  = mujacPart_qd + muMpart_qd' - muf'*bc_J*nu_qd;

  %wrt tau
  mujacPart_tau =modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,bc_tau_ext_MX);
  muMpart_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau);  %This is correct  
  ExtraParts_tau  = mujacPart_tau + muMpart_tau' - muf'*bc_J*nu_tau;

  %eqn 32 a 
  mu_rhs_q = jacobian([mu;muf]'*bc_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*bc_b,qd);
  
  %eqn 32 b 
  muHvec = [mu;muf]'*bc_K*aVec; 
  
  %Derivs 
  part31a_Brb = jacobian(mu_rhs_q,q); 
  part31b_Brb = jacobian(jacobian(muHvec,q),q);
  part31c_Brb = jacobian(ExtraParts_q,q); % jtimes([jac*nu_q]',q, muf);% + muf'*jac*nu_q,q); 
  bc_out_qq = part31a_Brb - (part31c_Brb+part31c_Brb') - part31b_Brb;
  
  %This is qqd  - works
  part31a_Brb = jacobian(mu_rhs_q,qd); 
  part31b_Brb = jacobian(jacobian(muHvec,q),qd); % + jacobian(jacobian(muf'*a_df_end,q),qd);
  part31c_Brb = jacobian(ExtraParts_qd,q); %+ muf'*jac*nu_q,q); 
  bc_out_qqd = part31a_Brb - (part31c_Brb)'-part31b_Brb;
  
  bc_out_qdqd = jacobian(mu_rhs_qd,qd); %Eqn 31a  
  bc_out_qtau = -jacobian(ExtraParts_tau,q);%+ muf'*jac*nu_tau,q);
  bc_secondMat  = [bc_out_qq bc_out_qqd; bc_out_qqd' bc_out_qdqd]; 
  
  FOP = [nu_q nu_qd nu_tau];
  f_Q = [f_q f_qd f_tau];
  bb=[mu;muf];
  bc_second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{bc_secondMat,bc_out_qtau}); 
  [bc_OutSecondMat, bc_OutQTau] = bc_second_mat_fn(q,qd,u,bc_dyn,bc_Kinv_vec_T(x,bb),FOP,f_Q); 
  bc_OutSecondMat_fn = Function('OutSecondMat_fn',{q,qd,u,bb,FOP,f_Q},{bc_OutSecondMat,bc_OutQTau}); 

%% Free Dynamics 
  %Eqn 32 part c
  %wrt q
  muMpart_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q);  %This is correct  
  ExtraParts_q  =  muMpart_q';

  %wrt qd
  muMpart_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd);  %This is correct  
  ExtraParts_qd  =  muMpart_qd';

  %wrt tau
  muMpart_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau);  %This is correct  
  ExtraParts_tau  =  muMpart_tau';

  %eqn 32 a 
  mu_rhs_q = jacobian(mu'*tau_b,q);
  mu_rhs_qd = jacobian(mu'*tau_b,qd);
  
  %eqn 32 b 
  muHvec = mu'*H*a7Vec; 
  
  %Derivs 
  part31a_Brb = jacobian(mu_rhs_q,q); 
  part31b_Brb = jacobian(jacobian(muHvec,q),q);
  part31c_Brb = jacobian(ExtraParts_q,q); % jtimes([jac*nu_q]',q, muf);% + muf'*jac*nu_q,q); 
  fr_out_qq = part31a_Brb - (part31c_Brb+part31c_Brb') - part31b_Brb;
  
  %This is qqd  - works
  part31a_Brb = jacobian(mu_rhs_q,qd); 
  part31b_Brb = jacobian(jacobian(muHvec,q),qd); % + jacobian(jacobian(muf'*a_df_end,q),qd);
  part31c_Brb = jacobian(ExtraParts_qd,q); %+ muf'*jac*nu_q,q); 
  fr_out_qqd = part31a_Brb - (part31c_Brb)'-part31b_Brb;
  
  fr_out_qdqd = jacobian(mu_rhs_qd,qd); %Eqn 31a  
  fr_out_qtau = -jacobian(ExtraParts_tau,q);%+ muf'*jac*nu_tau,q);
  fr_secondMat  = [fr_out_qq fr_out_qqd; fr_out_qqd' fr_out_qdqd]; 
  
  FOP = [nu_q nu_qd nu_tau];
  bb=mu;
  second_mat_fn = Function('second_mat_fn',{q,qd,u,a7Vec,bb,FOP},{fr_secondMat,fr_out_qtau}); 
  [fr_OutSecondMat, fr_OutQTau] = second_mat_fn(q,qd,u,fr_dyn,fr_Kinv_vec_T(x,bb),FOP); 
  fr_OutSecondMat_fn = Function('OutSecondMat_fn',{q,qd,u,bb,FOP},{fr_OutSecondMat,fr_OutQTau}); 

%% Front Impact Dynamics 
  %Eqn 32 part c
  %wrt q
  mujacPart_q = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_q_ext_MX);
  muMpart_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q);  %This is correct  
  ExtraParts_q  = mujacPart_q + muMpart_q' - muf'*ft_J*nu_q;

  %wrt qd
  mujacPart_qd = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_qd_ext_MX);
  muMpart_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd);  %This is correct  
  ExtraParts_qd  = mujacPart_qd + muMpart_qd' - muf'*ft_J*nu_qd;

  %wrt tau
  mujacPart_tau = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_tau_ext_MX);
  muMpart_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau);  %This is correct  
  ExtraParts_tau  = mujacPart_tau + muMpart_tau' - muf'*ft_J*nu_tau;

  %eqn 32 a 
  mu_rhs_q = jacobian([mu;muf]'*Im_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*Im_b,qd);
  
  %eqn 32 b 
  muHvec = [mu;muf]'*ft_K*aVec; 
  
  %Derivs 
  part31a_Brb = jacobian(mu_rhs_q,q); 
  part31b_Brb = jacobian(jacobian(muHvec,q),q);
  part31c_Brb = jacobian(ExtraParts_q,q); % jtimes([jac*nu_q]',q, muf);% + muf'*jac*nu_q,q); 
  ftIm_out_qq = part31a_Brb - (part31c_Brb+part31c_Brb') - part31b_Brb;
  
  %This is qqd  - works
  part31a_Brb = jacobian(mu_rhs_q,qd); 
  part31b_Brb = jacobian(jacobian(muHvec,q),qd); % + jacobian(jacobian(muf'*a_df_end,q),qd);
  part31c_Brb = jacobian(ExtraParts_qd,q); %+ muf'*jac*nu_q,q); 
  ftIm_out_qqd = part31a_Brb - (part31c_Brb)'-part31b_Brb;
  
  ftIm_out_qdqd = jacobian(mu_rhs_qd,qd); %Eqn 31a  
  ftIm_out_qtau = -jacobian(ExtraParts_tau,q);%+ muf'*jac*nu_tau,q);
  ftIm_secondMat  = [ftIm_out_qq ftIm_out_qqd; ftIm_out_qqd' ftIm_out_qdqd]; 
  
  FOP = [nu_q nu_qd nu_tau];
  f_Q = [f_q f_qd f_tau];
  bb=[mu;muf];
  second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{ftIm_secondMat,ftIm_out_qtau}); 
  [ftIm_OutSecondMat, ftIm_OutQTau] = second_mat_fn(q,qd,u,ft_Im_Dyn,ft_Kinv_vec_T(x,bb),FOP,f_Q); 
  ftIm_OutSecondMat_fn = Function('OutSecondMat_fn',{q,qd,u,bb,FOP,f_Q},{ftIm_OutSecondMat,ftIm_OutQTau}); 

%% Back Impact Dynamics 
  %Eqn 32 part c
  %wrt q
  %
  mujacPart_q = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,bc_q_ext_MX);
  muMpart_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q);  %This is correct  
  ExtraParts_q  = mujacPart_q + muMpart_q' - muf'*bc_J*nu_q;

  %wrt qd
  mujacPart_qd = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,bc_qd_ext_MX);
  muMpart_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd);  %This is correct  
  ExtraParts_qd  = mujacPart_qd + muMpart_qd' - muf'*bc_J*nu_qd;

  %wrt tau
  mujacPart_tau = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,bc_tau_ext_MX);
  muMpart_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau);  %This is correct  
  ExtraParts_tau  = mujacPart_tau + muMpart_tau' - muf'*bc_J*nu_tau;

  %eqn 32 a 
  %You can unify Back stance and Back impact setting mu_rhs_q as symbolic
  mu_rhs_q = jacobian([mu;muf]'*Im_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*Im_b,qd);
  
  %eqn 32 b 
  muHvec = [mu;muf]'*bc_K*aVec; 
  
  %Derivs 
  part31a_Brb = jacobian(mu_rhs_q,q); 
  part31b_Brb = jacobian(jacobian(muHvec,q),q);
  part31c_Brb = jacobian(ExtraParts_q,q); % jtimes([jac*nu_q]',q, muf);% + muf'*jac*nu_q,q); 
  bcIm_out_qq = part31a_Brb - (part31c_Brb+part31c_Brb') - part31b_Brb;
  
  %This is qqd  - works
  part31a_Brb = jacobian(mu_rhs_q,qd); 
  part31b_Brb = jacobian(jacobian(muHvec,q),qd); % + jacobian(jacobian(muf'*a_df_end,q),qd);
  part31c_Brb = jacobian(ExtraParts_qd,q); %+ muf'*jac*nu_q,q); 
  bcIm_out_qqd = part31a_Brb - (part31c_Brb)'-part31b_Brb;
  
  bcIm_out_qdqd = jacobian(mu_rhs_qd,qd); %Eqn 31a  
  bcIm_out_qtau = -jacobian(ExtraParts_tau,q);%+ muf'*jac*nu_tau,q);
  bcIm_secondMat  = [bcIm_out_qq bcIm_out_qqd; bcIm_out_qqd' bcIm_out_qdqd]; 
  
  FOP = [nu_q nu_qd nu_tau];
  f_Q = [f_q f_qd f_tau];
  bb=[mu;muf];
  second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{bcIm_secondMat,bcIm_out_qtau}); 
  %}
  [bcIm_OutSecondMat, bcIm_OutQTau] = second_mat_fn(q,qd,u,bc_Im_Dyn,bc_Kinv_vec_T(x,bb),FOP,f_Q); 
  bcIm_OutSecondMat_fn = Function('OutSecondMat_fn',{q,qd,u,bb,FOP,f_Q},{bcIm_OutSecondMat,bcIm_OutQTau}); 
   
%% EXTENDED MODIFIED RNEA
%%  Front Stance Dynamics 


    % come back:: Doing both at same time failed 
    mujacPart_qd = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_qd_ext_MX);
    muMpart_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd);        
    ExtraParts_qd  = mujacPart_qd + muMpart_qd';

    mujacPart_q = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_q_ext_MX);
    muMpart_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q);        
    ExtraParts_q  = mujacPart_q + muMpart_q';

    mujacPart_tau = modID_casadi_Extended(robot_no_grav,q,0*qd,0*qd,mu,ft_tau_ext_MX);
    muMpart_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau);        
    ExtraParts_tau  = mujacPart_tau + muMpart_tau';

    %Second partials of Extended Mod RNEA with q 
    A = jacobian(ExtraParts_q,q); %This is w in Eqn 32
    P =jacobian(-muf'*jac*nu_q,q); %it's really nu_q but see above comments:: This is u in Eqn 32
    down=hessian(-muf'*acc_act,q) + (P+P'); %Eqn 32a
    upper =  -hessian(out_Extended,q) - (A+A');%Eqn 32b
    ExtModRNEA_qq_V2 = upper + down;
    funcs.ExtModRNEA_qq_V2 = Function('ExtModRNEA_qq_V2',{q,qd,qdd,mu,muf,Fc_MX,f_q,nu_q},{ExtModRNEA_qq_V2});

    %Second partials of Extended Mod RNEA with qd
    upper = -hessian(out_Extended,qd) - jacobian(ExtraParts_qd,qd);%Eqn 32a
    down = hessian(-muf'*acc_act,qd) +jacobian(-muf'*jac*nu_qd,qd);%Eqn 32b
    ExtModRNEA_qdqd_V2 = upper + down; 
    funcs.ExtModRNEA_qdqd_V2= Function('ExtModRNEA_qdqd_V2',{q,qd,qdd,mu,muf,Fc_MX,f_qd,nu_qd},{ExtModRNEA_qdqd_V2});

    %Second partials of Extended Mod RNEA with q tau 
    upper =  -jacobian(ExtraParts_tau,q); %Eqn 32a 
    b1=jacobian(-muf'*jac*nu_tau,q); %Eqn 32b
    ExtModRNEA_qtau_V2 = upper + b1; 
    funcs.ExtModRNEA_qtau_V2 = Function('ExtModRNEA_qtau_V2',{q,mu,muf,Fc_MX,f_tau,nu_tau},{ExtModRNEA_qtau_V2});
    %here f_qd == f_tau and nu_qd = nu_tau

    %Second partials of Extended Mod RNEA with q qd
    upper  = -jacobian(jacobian(out_Extended,qd),q) - jacobian(ExtraParts_qd,q); %Eqn 32a 
    down=jacobian(jacobian(-muf'*acc_act,qd),q) + jacobian(-muf'*jac*nu_qd,q);%Eqn 32b 
    ExtModRNEA_qqd_V2 = upper + down; 
    funcs.ExtModRNEA_qqd_V2 = Function('ExtModRNEA_qqd_V2',{q,qd,qdd,mu,muf,Fc_MX,f_qd,nu_qd},{ExtModRNEA_qqd_V2}); 

    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2; 
            ExtModRNEA_qqd_V2' ExtModRNEA_qdqd_V2]; 

%}






   

  
 % You're HERE!!!!  Explicit Tensor stuff is done 
 % To do : via Mod RNEA 
 % To do : Tensor stuff 
  
  kkt_out = ft_dyn;
  
  Dyn_funcs.kkt_x = Function('KKT_x', {x,u},  { jacobian(jacobian(kkt_out, x),x) } );
  Dyn_funcs.kkt_qtau= Function('KKT_qtau', {q,qd,u},  { jacobian(jacobian(kkt_out, q),u) } );
  

  %Does it work 
  1==1;
  qVal = rand(size(q)); qdVal = rand(size(qd));
  tauVal = rand(size(u)); bbVal = rand(size(bb)); 
  xVal = [qVal; qdVal]; 
  h1 = Dyn_funcs.ft_dyn_x(xVal,tauVal); 
  h2 = Dyn_funcs.ft_dyn_u(xVal,tauVal);
  all = [h1 h2]; 
  FOPVal = all(1:Nb,:);
  fQVal = all(Nb+1:end,:);
        
    kkt_xx = full(Dyn_funcs.kkt_x(xVal,tauVal)); 
    kkt_xx = reshape(kkt_xx,[(Nb+2),Nb*2,Nb*2]); 
    v_xx_kkt = full(Tens3byVec(kkt_xx,bbVal','pre'));
   
    kkt_qtau = full(Dyn_funcs.kkt_qtau(qVal,qdVal,tauVal)); 
    kkt_qtau = reshape(kkt_qtau,[Nb+2 Nb length(u)]);
    v_qtau_kkt = (full(Tens3byVec(kkt_qtau,bbVal','pre'))); 
    
    bOther = v_qtau_kkt;
    [a b] = ft_OutSecondMat_fn(qVal,qdVal,tauVal,bbVal,FOPVal,fQVal);
    a = full(a - v_xx_kkt)
    b = full(b - bOther')
    
%     aVal = Dyn_funcs.ft_dyn(xVal,tauVal);
%     [GG, BB] = second_mat_fn(qVal,qdVal,tauVal,aVal,bbVal, FOPVal,fQVal)
    1==1;
    G
  



%{





if isa(q,'sym')
    ft_Jd = sym(zeros(size(ft_J)));             % jacobian derivates for two contact points
    bc_Jd = sym(zeros(size(bc_J)));
    for i = 1:size(ft_J,1)
        for j = 1:size(ft_J,2)
            ft_Jd(i,j) = jacobian(ft_J(i,j),q)*qd;
            bc_Jd(i,j) = jacobian(bc_J(i,j),q)*qd;
        end
    end
else 
   ft_Jd = jtimes(ft_J,q,qd); 
   bc_Jd = jtimes(bc_J,q,qd);
   
end

% front stance dynamics
ft_b = [tau_b;ft_Jd*qd];

%%
if ~isa(q,'sym')
    %Acceleration of contact point
    [~,acc,Xup, v]= ID_casadi_Extended(robot,q,qd,0*qd);        
    %Acceleration of contact point expressed in inertia frame
    [jacdot_qd,~]= Fncontact_acc(acc,v,Xup,hh,[5 7],1);
    %Remove the y component 
    ft_jacdot_qd = [jacdot_qd{1}(1) jacdot_qd{1}(3)];
    bck_jacdot_qd = [jacdot_qd{2}(1) jacdot_qd{2}(3)];
    %This doesn't work yet! Ctrl+F "Tag"
    Dyn_funcs.bck_jacdot_qd = Function('bck_jacdot_qd',{x},{bck_jacdot_qd});
    Dyn_funcs.ft_jacdot_qd = Function('ft_jacdot_qd',{x},{ft_jacdot_qd}); 
end
%%
ft_b_x = jacobian(ft_b,x);
ft_b_u = jacobian(ft_b,u);
ft_K = [H, ft_J';
        -ft_J, zeros(size(ft_J,1),size(ft_J,1))]; % KKT Dynamics matrix
if isa(q,'sym')
    ft_K_x = MatrixPartialD(ft_K, x);         % size 9 by 9 by 14 tensor
else 
    ft_K_x = jacobian(ft_K,x);%needs to be reshaped
end 


% back stance dynamics
bc_b = [tau_b;bc_Jd*qd];
bc_b_x = jacobian(bc_b,x);
bc_b_u = jacobian(bc_b,u);
bc_K = [H, bc_J';
        -bc_J, zeros(size(bc_J,1),size(bc_J,1))]; % KKT Dynamics matrix
if isa(q,'sym')
    bc_K_x = MatrixPartialD(bc_K, x);         % size 9 by 9 by 14 tensor
else
    bc_K_x = jacobian(bc_K,x);
end

qd_x =  jacobian(qd,x);

% Impact Dynamics
a = [H*qd; zeros(size(ft_Jd*qd))]; % x-
a_x = jacobian(a,x);
a_u = jacobian(a,u); % all zeros
q_x = jacobian(q,x); % unify identity matrix and zeros
1==1;
if isa(q,'sym')
% matlabFunction(H, tau_b, H_x, tau_bx, tau_bu, qd_x, 'file', 'support/FreeDynamics', 'vars', {x, u});
% matlabFunction(ft_K,ft_b,ft_K_x,ft_b_x, ft_b_u, qd_x, 'file','support/FrontStanceDyn','vars',{x, u});
% matlabFunction(bc_K,bc_b,bc_K_x,bc_b_x, bc_b_u, qd_x, 'file','support/BackStanceDyn','vars',{x, u});
% 
% matlabFunction(ft_K,a,ft_K_x,a_x, a_u, q_x, 'file','support/FrontImpactDyn','vars',{x, u});
% matlabFunction(bc_K,a,bc_K_x,a_x, a_u, q_x, 'file','support/BackImpactDyn','vars',{x, u});
% matlabFunction(ft_J, ft_Jd,'file','support/FrontJacobian','vars',{x});
% matlabFunction(bc_J, bc_Jd,'file','support/BackJacobian','vars',{x});
else 
    Dyn_funcs.FreeDyn = Function('FreeDyn',{x,u},{H, tau_b, H_x, tau_bx, tau_bu, qd_x});
    Dyn_funcs.FrntStncDyn = Function('FrntStncDyn',{x,u},{ft_K,ft_b,ft_K_x,ft_b_x, ft_b_u, qd_x});
    Dyn_funcs.BckStncDyn = Function('BckStncDyn',{x,u},{bc_K,bc_b,bc_K_x,bc_b_x, bc_b_u, qd_x});
    
    Dyn_funcs.FrntImpctDyn = Function('FrntImpctDyn',{x,u},{ft_K,a,ft_K_x,a_x, a_u, q_x}); 
    Dyn_funcs.BckImpctDyn =  Function('BckImpctDyn',{x,u},{bc_K,a,bc_K_x,a_x, a_u, q_x});
    Dyn_funcs.FrntJac = Function('FrntJac',{x},{ft_J, ft_Jd}); 
    Dyn_funcs.BckJac = Function('BckJac',{x},{bc_J, bc_Jd});
end
end

function [a_cnt,v_cnt]= Fncontact_acc(a_jnt,v_jnt,Xup,pos,stateIndx, gravComp) 
%returns acceleration at contact points
%arguments:: (a_jnt,v_jnt,Xup,Nb,gravComp) 
%Returns:: a_contact = jacdot*qdd and v_contact
sz = length(stateIndx); 
a_cnt = cell(1,sz); v_cnt = a_cnt; acc_coriolis = a_cnt;

for idx = 1:sz
    a_cnt{idx} = Vpt(a_jnt{stateIndx(idx)},pos(:,idx)); %Return ax ay at contact
    v_cnt{idx} = Vpt(v_jnt{idx},pos(:,idx));

    %Corriolis effects of gravity
    acc_coriolis{idx} = cross( [0 0 v_jnt{stateIndx(idx)}(end)]',[v_cnt{idx}]);
end
X=eye(3); Xd = cell(size(Xup));
for i=1:length(Xup)
   X= Xup{i}(1:3,1:3)*X; %Rot to ineratial frame 
   Xd{i} = X;
end
for idx =1:sz 
   X = plnr(0,pos(:,idx))* Xd{stateIndx(idx)};
   [R_end, ~] = plux_plnr_v2(X); % Need modified version of plux that supports planar
                      % Compare position kinematics
   a_cnt{idx} = R_end'*(a_cnt{idx}+acc_coriolis{idx}); 
   if gravComp
        a_cnt{idx}  = a_cnt{idx} + [0 0 -9.81]'; 
   end
   v_cnt{idx} = R_end'*v_cnt{idx}; % Similar thing for velocities.                   
end

%}
end

