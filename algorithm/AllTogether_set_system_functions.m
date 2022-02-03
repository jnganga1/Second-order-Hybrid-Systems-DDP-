function [params, robot_params,Dyn_funcs] = AllTogether_set_system_functions(x, u, mode)
addpath([pwd '/support']);
addpath([pwd '/algorithm']);
import casadi.*

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
robot_params.g = 1.62; 

[H, C] = HandC_casadi(robot, q, qd);

robot_params.x_size = x_size;
robot_params.u_size = u_size;
robot_params.q_size = q_size;
robot_params.y_size = y_size;

%% Dynamics Formulation
% Free Dynamics (flying-phase dynamics)
S = [zeros(q_size-u_size, u_size);
     eye(u_size)]; % selection matrix
tau_b = S*u - C;    % control term, gravity term, centrifugal term and corriolis term

% Contact jacobian and its derivative
ft = 1; bc = 2;
Pos_ctact = fkin_quadruped_2D(q, robot_params);  % forward kinematics (two contact points)
hh = Pos_ctact;
Pos_ctact = [Pos_ctact(1,:); Pos_ctact(3,:)];% eliminate y (2D quadruped in x-z plane)

%Our Rotations happen in 3d space, adding all zero y position
Pos_ctact_3d = [hh(1,:); 0*hh(2,:) ; hh(3,:)]; 

ft_J = jacobian(Pos_ctact(:,ft),q);
bc_J = jacobian(Pos_ctact(:,bc),q);

%Time derivatives of the jacobian
ft_Jd = jtimes(ft_J,q,qd); 
bc_Jd = jtimes(bc_J,q,qd);

Dyn_funcs.FrntJac = Function('FrntJac',{x},{ft_J, ft_Jd}); 
Dyn_funcs.BckJac = Function('BckJac',{x},{bc_J, bc_Jd});

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


%% Symbolics and Functions for later

%Can replace the right hand of any mode, i.e., aVec = K^(-1) * rhs
aVec = MX.sym('aVec',[9,1]);
%Can replace the right hand of the free mode, i.e., aVec = H^(-1) * rhs
a7Vec = MX.sym('aVec',[7,1]); 

%For simplicity, redefine
Fc_MX  = aVec(Nb+1:end);
qdd = aVec(1:Nb,:);
%In case you need a virtuous force in 3d 
Fc_MX_3d  = [Fc_MX(1,:);MX.zeros(1,1);Fc_MX(2,:)];

%Will become the qqd part of gamma^T * ()
mu = MX.sym('mu',Nb); 
%Will become the force part of gamma^T * ()
muf = MX.sym('muf',[2 1]); 
      
%The qdd parts of first_qMX aka the first derivatives of q,qd,tau
nu_q = MX.sym('nu_q',[Nb,Nb]);
nu_qd = MX.sym('nu_qd',[Nb,Nb]);
nu_tau = MX.sym('nu_tau',[Nb,u_size]);

%The lambda parts of first_qMX aka the first derivatives of q,qd,tau
f_q  = MX.sym('f_q',[2,Nb]);
f_qd  = MX.sym('f_qd',[2,Nb]); 
f_tau = MX.sym('f_tau',[2,u_size]); 

%Zeroing out the y dimension since we are in the 2d plane 
f_Q  = [f_q(1,:);MX.zeros(1,7);f_q(2,:)];
f_Qd  = [f_qd(1,:);MX.zeros(1,7);f_qd(2,:)];
f_Tau = [f_tau(1,:);MX.zeros(1,u_size);f_tau(2,:)];

%Initialize with the proper size
ft_spat_f_q = MX.zeros(6,7); bc_spat_f_q = MX.zeros(6,7);
ft_spat_f_qd = MX.zeros(6,7); bc_spat_f_qd= MX.zeros(6,7);
ft_spat_f_tau = MX.zeros(6,u_size); bc_spat_f_tau = MX.zeros(6,u_size);
for idx = 1:Nb 
    %Fpt assumes a (x,y,z) point and an (fx,fy,fz) force. Since our "Force"
    %is a matrix (i.e. partials), we recursively use the columns as forces.
    ft_spat_f_q(:,idx) = Fpt(-f_Q(:,idx),Pos_ctact_3d(:,ft)); %Introduces jacobian to f_q  
    ft_spat_f_qd(:,idx) = Fpt(-f_Qd(:,idx),Pos_ctact_3d(:,ft)); %Introduces jacobian to f_q 
    
    bc_spat_f_q(:,idx) = Fpt(-f_Q(:,idx),Pos_ctact_3d(:,bc)); %Introduces jacobian to f_q  
    bc_spat_f_qd(:,idx) = Fpt(-f_Qd(:,idx),Pos_ctact_3d(:,bc)); %Introduces jacobian to f_q 
end
for idx = 1:u_size
    ft_spat_f_tau(:,idx) = Fpt(-f_Tau(:,idx),Pos_ctact_3d(:,ft)); %Introduces jacobian to f_q  
    bc_spat_f_tau(:,idx) = Fpt(-f_Tau(:,idx),Pos_ctact_3d(:,bc)); %Introduces jacobian to f_q  
end 

%Deriv in format Required by RNEA algo ::state 5 is front legs. state 7 is hind legs
ft_q_ext_MX = cell([1 Nb]); ft_q_ext_MX{5} = ft_spat_f_q;
ft_qd_ext_MX = cell([1 Nb]); ft_qd_ext_MX{5} = ft_spat_f_qd;
ft_tau_ext_MX = cell([1 Nb]); ft_tau_ext_MX{5} = ft_spat_f_tau;

bc_q_ext_MX = cell([1 Nb]); bc_q_ext_MX{7} = bc_spat_f_q;
bc_qd_ext_MX = cell([1 Nb]); bc_qd_ext_MX{7} = bc_spat_f_qd;
bc_tau_ext_MX = cell([1 Nb]); bc_tau_ext_MX{7} = bc_spat_f_tau;
 
% Forces in format Required by RNEA algo ::state 5 is front legs. state 7 is hind legs
ft_force_q__MX = cell([1 Nb]); ft_force_q__MX{5} = Fpt(-Fc_MX_3d,Pos_ctact_3d(:,ft));
bc_force_q__MX = cell([1 Nb]); bc_force_q__MX{7} = Fpt(-Fc_MX_3d,Pos_ctact_3d(:,bc));

%zero out the robot's gravity 
robot_no_grav = robot; robot_no_grav.gravity = 0*robot_no_grav.gravity;

%% Functions
%This functions does ft_K * aVec where aVec is any 9x1 vector

[A,a_jnt,~,v_jnt]=ID_casadi_Extended(robot_no_grav,q,qd*0,qdd,ft_force_q__MX);
ftJqdd = jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
ft_K_vec = Function('ft_K_vec',{q,aVec},...
    {[A; -ftJqdd]}); 
%This functions does bc_K * aVec where aVec is any 9x1 vector
[A,a_jnt,~,v_jnt]=ID_casadi_Extended(robot_no_grav,q,qd*0,qdd,bc_force_q__MX);
bcJqdd = jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
bc_K_vec = Function('bc_K_vec',{q,aVec},...
    {[A ; -bcJqdd]}); 
%This functions does H * a7Vec where a7Vec is any 7x1 vector
H_vec = Function('ft_K_vec',{q,a7Vec},...
    {ID_casadi(robot_no_grav,q,qd*0,a7Vec)});


%%
%% This is a place holder for the later Hinv Function -> This will soon be replaced with a function call 

%Explicit formation of Free Hinv Matrix
% Hinv = H \ a7Vec; 
% fr_Kinv_vec = Function('fr_Hinv_vec',{x,a7Vec},{Hinv});
% ABA to do H inverse 
fr_Kinv_vec = Function('fr_Kinv_vec',{x,a7Vec}, {FDab_casadi(robot_no_grav, q, qd*0, a7Vec)}); 

Hinv_T = a7Vec' / H; 
fr_Kinv_vec_T = Function('fr_Hinv_vec',{x,a7Vec},{Hinv_T});

ft_Kinv = ft_K \ aVec; 
ft_Kinv_vec = Function('ft_Kinv_vec',{x,aVec},{ft_Kinv}); 

ft_Kinv_T =  aVec' / ft_K; 
ft_Kinv_vec_T = Function('ft_Kinv_vec',{x,aVec},{ft_Kinv_T}); 

% Out = Hinv_Vec(x,robot,robot_params,ft_force_q__MX,bc_force_q__MX);
% xVal = rand(size(x)); yVal =rand(size(aVec));
% a = full(Out(xVal,yVal))'
% b = full(ft_Kinv_vec_T(xVal,yVal))
% a-b


bc_Kinv = bc_K \ aVec; 
bc_Kinv_vec = Function('bc_Kinv_vec',{x,aVec},{bc_Kinv });

bc_Kinv_T =  aVec' / bc_K; 
bc_Kinv_vec_T = Function('bc_Kinv_vec',{x,aVec},{bc_Kinv_T});

%%

%Front Stance Dynamics
% ft_dyn = ft_K\ft_b; 
ft_dyn = ft_Kinv_vec(x,ft_b); 
Dyn_funcs.ft_dyn = Function('ft_dyn',{x,u},{ft_dyn});
ft_dyn_x_w_avec = ft_Kinv_vec(x,jacobian(ft_b,x) - jacobian(ft_K_vec(q,aVec),x));
ft_dyn_x_w_avec_fn = Function('ft_dyn_x_w_avec_fn',{x,aVec},{ft_dyn_x_w_avec});
ft_dyn_x_soln = ft_dyn_x_w_avec_fn(x,ft_dyn);
ft_dyn_u = ft_Kinv_vec(x,jacobian(ft_b,u));
% Dyn_funcs.ft_dyn_u = Function('ft_dyn_u',{x,u},{ft_dyn_u});
all_ft_dyn = [ft_dyn_x_soln,ft_dyn_u]; 
Dyn_funcs.ft_dyn_derivs = Function('ft_dyn_derivs',{x,u},{all_ft_dyn});
Dyn_funcs.Ignore = Function('ag',{x,u},{ft_dyn_x_soln});
FOP_ft_dyn = all_ft_dyn(1:q_size,:); 
FQ_ft_dyn = all_ft_dyn(q_size+1:end,:);

%Back Stance Dynamics
% bc_dyn = bc_K \bc_b; 
bc_dyn = bc_Kinv_vec(x,bc_b); 
Dyn_funcs.bc_dyn = Function('bc_dyn',{x,u},{bc_dyn});
bc_dyn_x_w_avec = bc_Kinv_vec(x,jacobian(bc_b,x) - jacobian(bc_K_vec(q,aVec),x));
bc_dyn_x_w_avec_fn = Function('bc_dyn_x_w_avec_fn',{x,aVec},{bc_dyn_x_w_avec});
bc_dyn_x_soln = bc_dyn_x_w_avec_fn(x,bc_dyn); 
bc_dyn_u = bc_Kinv_vec(x,jacobian(bc_b,u));
all_bc_dyn = [bc_dyn_x_soln,bc_dyn_u]; 
Dyn_funcs.bc_dyn_derivs = Function('bc_dyn_derivs',{x,u},{all_bc_dyn});
FOP_bc_dyn = all_bc_dyn(1:q_size,:); 
FQ_bc_dyn = all_bc_dyn(q_size+1:end,:);

%Free Dynamics
fr_dyn = fr_Kinv_vec(x,tau_b); 
Dyn_funcs.fr_dyn = Function('fr_dyn',{x,u},{fr_dyn});
fr_dyn_x_w_avec = fr_Kinv_vec(x,jacobian(tau_b,x) - jacobian(H_vec(q,a7Vec),x));
fr_dyn_x_w_avec_fn = Function('fr_dyn_x_w_avec_fn',{x,a7Vec},{fr_dyn_x_w_avec});
fr_dyn_x_soln = fr_dyn_x_w_avec_fn(x,fr_dyn); 
fr_dyn_u = fr_Kinv_vec(x,jacobian(tau_b,u));
all_fr_dyn = [fr_dyn_x_soln,fr_dyn_u];
Dyn_funcs.fr_dyn_derivs = Function('fr_dyn_derivs',{x,u},{all_fr_dyn});
FOP_fr_dyn = all_fr_dyn; %just naming convention. ignore

%Front Impact Dynamics h
% ft_Im_Dyn = ft_K\Im_b; 
ft_Im_Dyn = ft_Kinv_vec(x,Im_b);
Dyn_funcs.ft_Im_Dyn = Function('fr_Im_Dyn',{x},{ft_Im_Dyn});
ftIm_dyn_x_w_avec = ft_Kinv_vec(x,jacobian(Im_b,x) - jacobian(ft_K_vec(q,aVec),x));
ftIm_dyn_x_w_avec_fn = Function('ftIm_dyn_x_w_avec_fn',{x,aVec},{ftIm_dyn_x_w_avec});
ftIm_dyn_x_soln = ftIm_dyn_x_w_avec_fn(x,ft_Im_Dyn); 
ftIm_dyn_u = ft_Kinv_vec(x,jacobian(Im_b,u));
all_ft_Im_dyn = [ftIm_dyn_x_soln,ftIm_dyn_u]; 
Dyn_funcs.ft_Im_dyn_derivs = Function('ft_dyn_u',{x,u},{all_ft_Im_dyn});
FOP_ft_Im_dyn = all_ft_Im_dyn(1:q_size,:); 
FQ_ft_Im_dyn = all_ft_Im_dyn(q_size+1:end,:);


%Back Impact Dynamics 
% bc_Im_dyn = bc_K \Im_b; 
bc_Im_dyn  = bc_Kinv_vec(x,Im_b);
Dyn_funcs.bc_Im_dyn = Function('bc_Im_dyn',{x},{bc_Im_dyn});
bc_Im_dyn_x_w_avec = bc_Kinv_vec(x,jacobian(Im_b,x) - jacobian(bc_K_vec(q,aVec),x));
bc_Im_dyn_x_w_avec_fn = Function('ftIm_dyn_x_w_avec_fn',{x,aVec},{bc_Im_dyn_x_w_avec});
bc_Im_dyn_x_soln = bc_Im_dyn_x_w_avec_fn(x,bc_Im_dyn); 
bc_Im_dyn_u = bc_Kinv_vec(x,jacobian(Im_b,u)); %This is all zeros 
all_bc_Im_dyn = [bc_Im_dyn_x_soln,bc_Im_dyn_u];
Dyn_funcs.bc_Im_dyn_derivs = Function('ft_dyn_u',{x,u},{all_bc_Im_dyn});
FOP_bc_Im_dyn = all_bc_Im_dyn(1:q_size,:); 
FQ_bc_Im_dyn = all_bc_Im_dyn(q_size+1:end,:);

%Verification Process 
%  %add { here to comment out 
Dyn_funcs.kkt_ft_first = Function('ft_first',{x,u},{
    [jacobian(ft_dyn,x),jacobian(ft_dyn,u)]});
Dyn_funcs.kkt_bc_first = Function('bc_first',{x,u},{
    [jacobian(bc_dyn,x),jacobian(bc_dyn,u)]});
Dyn_funcs.kkt_fr_first = Function('fr_first',{x,u},{
    [jacobian(fr_dyn,x),jacobian(fr_dyn,u)]});
Dyn_funcs.kkt_ft_Im_first = Function('ft_Im_first',{x,u},{
    [jacobian(ft_Im_Dyn,x),jacobian(ft_Im_Dyn,u)]});
Dyn_funcs.kkt_bc_Im_first = Function('bc_Im_first',{x,u},{
    [jacobian(bc_Im_dyn,x),jacobian(bc_Im_dyn,u)]});


%}

%%
%% Second order 
%% 
% Might reformat this so that each method is its own script -> would be
% cleaner
% if SecondOrder
%% EXPLICIT METHODS
%% Front Stance Dynamics 
  %Eqn 32 part c
  %wrt q
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,ft_q_ext_MX); 
  ftJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
   ExtraParts_q = A - ftJnuq;
   
  %wrt qd
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,ft_qd_ext_MX); 
  ftJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
  ExtraParts_qd = A - ftJnuqd;

  %wrt tau
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,ft_tau_ext_MX); 
  ftJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
  ExtraParts_tau = A - ftJnutau;

  %eqn 32 a 
  mu_rhs_q = jacobian([mu;muf]'*ft_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*ft_b,qd);
  
  %eqn 32 b 
%   muHvec = [mu;muf]'*ft_K*aVec; 
  % Returns mu * [H qdd + ft_J *f]
 [muHvec_v2,a_jnt,~,v_jnt] =modID_casadi_Extended(robot_no_grav,q,qd*0,qdd,mu,ft_force_q__MX);
 [ftqdd] = -jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
  muHvec = muHvec_v2 + muf' * ftqdd; % This is muf * -ft_J * qdd;

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
  [ft_OutSecondMat, ft_OutQTau] = second_mat_fn(q,qd,u,ft_dyn,ft_Kinv_vec_T(x,bb),FOP_ft_dyn,FQ_ft_dyn); 
  ft_OutSecondMat_fn = Function('OutSecondMat_fn',{x,u,bb},{all_ft_dyn,ft_OutSecondMat,ft_OutQTau}); 

%% Back Stance Dynamics 
  %Eqn 32 part c
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,bc_q_ext_MX); %H*nu_q + bc_J'*f_q
  bcJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0); %bc_J'*nu_q
  ExtraParts_q = A - bcJnuq;
  
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,bc_qd_ext_MX); 
  bcJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
  ExtraParts_qd = A - bcJnuqd;
  
  %wrt tau 
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,bc_tau_ext_MX); 
  bcJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
  ExtraParts_tau = A - bcJnutau;

  %eqn 32 a 
  mu_rhs_q = jacobian([mu;muf]'*bc_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*bc_b,qd);
  
  %eqn 32 b 
%   muHvec = [mu;muf]'*bc_K*aVec; 
  % Returns mu * [H qdd + bc_K *f]
  [muHvec_v2,a_jnt,~,v_jnt] =modID_casadi_Extended(robot_no_grav,q,qd*0,qdd,mu,bc_force_q__MX);
 [bcqdd] = -jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
  muHvec = muHvec_v2 + muf' * bcqdd;  
  
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
  second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{bc_secondMat,bc_out_qtau});
  
  [bc_OutSecondMat, bc_OutQTau] = second_mat_fn(q,qd,u,bc_dyn,bc_Kinv_vec_T(x,bb),FOP_bc_dyn,FQ_bc_dyn); 
  bc_OutSecondMat_fn = Function('OutSecondMat_fn',{x,u,bb},{all_bc_dyn,bc_OutSecondMat,bc_OutQTau}); 

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
%   muHvec = mu'*H*a7Vec; 
  muHvec = modID_casadi(robot_no_grav,q,0*qd,mu,a7Vec);  %This is correct  

  
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
  [fr_OutSecondMat, fr_OutQTau] = second_mat_fn(q,qd,u,fr_dyn,fr_Kinv_vec_T(x,bb),FOP_fr_dyn); 
  fr_OutSecondMat_fn = Function('OutSecondMat_fn',{x,u,bb},{all_fr_dyn,fr_OutSecondMat,fr_OutQTau}); 

%% Front Impact Dynamics 
  %Eqn 32 part c
  %wrt q
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,ft_q_ext_MX); 
  ftJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
  ExtraParts_q = A - ftJnuq;
 
  %wrt qd
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,ft_qd_ext_MX); 
  ftJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
  ExtraParts_qd = A - ftJnuqd;

  %wrt tau
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,ft_tau_ext_MX); 
  ftJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
  ExtraParts_tau= A - ftJnutau;
  
  %eqn 32 a 
  mu_rhs_q = jacobian([mu;muf]'*Im_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*Im_b,qd);
  
  %eqn 32 b 
%   muHvec = [mu;muf]'*ft_K*aVec; 
  % Returns mu * [H qdd + ft_J *f]
  [muHvec_v2,a_jnt,~,v_jnt] =modID_casadi_Extended(robot_no_grav,q,qd*0,qdd,mu,ft_force_q__MX);
  [ftqdd] = -jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
  muHvec = muHvec_v2 + muf' * ftqdd; 
  
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
  [ftIm_OutSecondMat, ftIm_OutQTau] = second_mat_fn(q,qd,u,ft_Im_Dyn,ft_Kinv_vec_T(x,bb),FOP_ft_Im_dyn,FQ_ft_Im_dyn); 
  ft_Im_OutSecondMat_fn = Function('OutSecondMat_fn',{x,u,bb},{all_ft_Im_dyn,ftIm_OutSecondMat,ftIm_OutQTau}); 

%% Back Impact Dynamics 
  %Eqn 32 part c
  %wrt q
  %
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,bc_q_ext_MX); 
  bcJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
  ExtraParts_q= A - bcJnuq;  

  %wrt qd
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,bc_qd_ext_MX); 
  bcJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
  ExtraParts_qd = A - bcJnuqd;

  %wrt tau 
  [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,bc_tau_ext_MX); 
  bcJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
  ExtraParts_tau= A - bcJnutau;

  %eqn 32 a 
  %You can unify Back stance and Back impact setting mu_rhs_q as symbolic
  mu_rhs_q = jacobian([mu;muf]'*Im_b,q);
  mu_rhs_qd = jacobian([mu;muf]'*Im_b,qd);
  
  %eqn 32 b 
%   muHvec = [mu;muf]'*bc_K*aVec; 
  % Returns mu * [H qdd + bc_K *f]
  [muHvec_v2,a_jnt,~,v_jnt] =modID_casadi_Extended(robot_no_grav,q,qd*0,qdd,mu,bc_force_q__MX);
  [bcqdd] = -jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
  muHvec = muHvec_v2 + muf' * bcqdd; 
  
  %Derivs 
  part31a_Brb = jacobian(mu_rhs_q,q); 
  part31b_Brb = jacobian(jacobian(muHvec,q),q);
  part31c_Brb = jacobian(ExtraParts_q,q); % jtimes([jac*nu_q]',q, muf);% + muf'*jac*nu_q,q); 
  bc_Im_out_qq = part31a_Brb - (part31c_Brb+part31c_Brb') - part31b_Brb;
  
  %This is qqd  - works
  part31a_Brb = jacobian(mu_rhs_q,qd); 
  part31b_Brb = jacobian(jacobian(muHvec,q),qd); % + jacobian(jacobian(muf'*a_df_end,q),qd);
  part31c_Brb = jacobian(ExtraParts_qd,q); %+ muf'*jac*nu_q,q); 
  bc_Im_out_qqd = part31a_Brb - (part31c_Brb)'-part31b_Brb;
  
  bc_Im_out_qdqd = jacobian(mu_rhs_qd,qd); %Eqn 31a  
  bc_Im_out_qtau = -jacobian(ExtraParts_tau,q);%+ muf'*jac*nu_tau,q);
  bc_Im_secondMat  = [bc_Im_out_qq bc_Im_out_qqd; bc_Im_out_qqd' bc_Im_out_qdqd]; 
  
  FOP = [nu_q nu_qd nu_tau];
  f_Q = [f_q f_qd f_tau];
  bb=[mu;muf];
  second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{bc_Im_secondMat,bc_Im_out_qtau}); 
  
  [bc_Im_OutSecondMat, bc_Im_OutQTau] = second_mat_fn(q,qd,u,bc_Im_dyn,bc_Kinv_vec_T(x,bb),FOP_bc_Im_dyn,FQ_bc_Im_dyn); 
  bc_Im_OutSecondMat_fn = Function('OutSecondMat_fn',{x,u,bb},{all_bc_Im_dyn,bc_Im_OutSecondMat,bc_Im_OutQTau}); 
  %} 
%% EXTENDED MODIFIED RNEA
%%  Front Stance Dynamics 
    [out_Extended,a_jnt,~,v_jnt] =modID_casadi_Extended(robot,q,qd,qdd,mu,ft_force_q__MX);
%     acc_act = (ft_J*qdd + ft_Jd * qd);  
    acc_act = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,1);

    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,ft_qd_ext_MX); 
    ftJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
    ExtraParts_qd = A - ftJnuqd;

    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,ft_q_ext_MX); 
    ftJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
    ExtraParts_q = A - ftJnuq;    
    
    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,ft_tau_ext_MX); 
    ftJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
    ExtraParts_tau = A - ftJnutau;
    
    %Second partials of Extended Mod RNEA with q
    d=hessian(-out_Extended + acc_act,q);
    b=jacobian(-ExtraParts_q,q);
    ExtModRNEA_qq_V2 = d + b + b';
    
    %Second partials of Extended Mod RNEA with qd
    upper = hessian(out_Extended -acc_act,qd) + jacobian(ExtraParts_qd,qd);%Eqn 32a    
    ExtModRNEA_qdqd_V2 = -upper;

    %Second partials of Extended Mod RNEA with q tau 
    upper =  -jacobian(ExtraParts_tau,q); %Eqn 32a 
    ExtModRNEA_qtau_V2 = upper;
  
    %Second partials of Extended Mod RNEA with q qd
    upper  = jacobian(jacobian(-out_Extended+acc_act,qd),q) + jacobian(-ExtraParts_qd,q); %Eqn 32a 
    ExtModRNEA_qqd_V2 = upper; 
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{SecondMat,ExtModRNEA_qtau_V2}); 
    [ft_OutSecondMat, ft_OutQTau] = second_mat_fn(q,qd,u,ft_dyn,ft_Kinv_vec_T(x,bb),FOP_ft_dyn,FQ_ft_dyn); 
    ft_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,bb},{all_ft_dyn,ft_OutSecondMat,ft_OutQTau}); 
%%  Back Stance Dynamics 
    [out_Extended,a_jnt,~,v_jnt] =modID_casadi_Extended(robot,q,qd,qdd,mu,bc_force_q__MX);
    acc_act = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,1);
    
   [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,bc_qd_ext_MX); 
    bcJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
    ExtraParts_qd = A - bcJnuqd;  
    
    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,bc_q_ext_MX); 
    bcJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
    ExtraParts_q = A - bcJnuq;
    
    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,bc_tau_ext_MX); 
    bcJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
    ExtraParts_tau= A - bcJnutau;

    %Second partials of Extended Mod RNEA with q
    d=hessian(-out_Extended + acc_act,q);
    b=jacobian(-ExtraParts_q ,q);
    ExtModRNEA_qq_V2 = d + b + b';
    
    %Second partials of Extended Mod RNEA with qd
    upper = hessian(out_Extended-acc_act,qd) + jacobian(ExtraParts_qd,qd);%Eqn 32a    
    ExtModRNEA_qdqd_V2 = -upper;

    %Second partials of Extended Mod RNEA with q tau 
    upper =  -jacobian(ExtraParts_tau,q); %Eqn 32a 
    ExtModRNEA_qtau_V2 = upper;
  
    %Second partials of Extended Mod RNEA with q qd
    upper  = jacobian(jacobian(-out_Extended+acc_act,qd),q) + jacobian(-ExtraParts_qd,q); %Eqn 32a 
    ExtModRNEA_qqd_V2 = upper; 
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{SecondMat,ExtModRNEA_qtau_V2}); 
    [bc_OutSecondMat, bc_OutQTau] = second_mat_fn(q,qd,u,bc_dyn,bc_Kinv_vec_T(x,bb),FOP_bc_dyn,FQ_bc_dyn); 
    bc_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,bb},{all_bc_dyn,bc_OutSecondMat,bc_OutQTau}); 
%%  Free Dynamics
%
    out_Extended =modID_casadi(robot,q,qd,a7Vec,mu);

    muMpart_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd);        
    ExtraParts_qd  = muMpart_qd';

    muMpart_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q);        
    ExtraParts_q  =  muMpart_q';

    muMpart_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau);        
    ExtraParts_tau  = muMpart_tau';

    %Second partials of Extended Mod RNEA with q
    d=hessian(-out_Extended,q);
    b=jacobian(-ExtraParts_q ,q);
    ExtModRNEA_qq_V2 = d + b + b';

    %Second partials of Extended Mod RNEA with qd
    upper = hessian(out_Extended,qd) + jacobian(ExtraParts_qd,qd);%Eqn 32a    
    ExtModRNEA_qdqd_V2 = -upper;

    %Second partials of Extended Mod RNEA with q tau 
    upper =  -jacobian(ExtraParts_tau,q); %Eqn 32a 
    ExtModRNEA_qtau_V2 = upper;

    %Second partials of Extended Mod RNEA with q qd
    upper  = jacobian(jacobian(-out_Extended,qd),q) + jacobian(-ExtraParts_qd,q); %Eqn 32a 
    ExtModRNEA_qqd_V2 = upper; 
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{x,u,a7Vec,mu,FOP},{SecondMat,ExtModRNEA_qtau_V2}); 
    [fr_OutSecondMat, fr_OutQTau] = second_mat_fn(x,u,fr_dyn,fr_Kinv_vec_T(x,mu),FOP_fr_dyn); 
    fr_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,mu},{all_fr_dyn,fr_OutSecondMat,fr_OutQTau}); 

%% Back Impact Dynamics
%
    [out_Extended,a_jnt,~,v_jnt] =modID_casadi_Extended(robot_no_grav,q,qd*0,qdd,mu,bc_force_q__MX); %H*qdd - bc_j*Fc
    out_Extended = -out_Extended + modID_casadi(robot_no_grav,q,0*qd,mu,qd); % mu'*H*qd;
    
    acc_act = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);    
    
    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,bc_q_ext_MX); 
    bcJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
    ExtraParts_q = A - bcJnuq;
    
    dd = jacobian(ExtraParts_q ,q);
    upper = hessian(out_Extended + acc_act,q) - (dd + dd');%Eqn 32a
    
    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,bc_qd_ext_MX); 
    bcJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
    ExtraParts_qd = A - bcJnuqd;

     [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,bc_tau_ext_MX); 
    bcJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,bc,0);
    ExtraParts_tau = A - bcJnutau;
       
    %Second partials of Extended Mod RNEA with qq
    ExtModRNEA_qq_V2 =upper;
    
    %Second partials of Extended Mod RNEA with qd
    upper = hessian(out_Extended+acc_act,qd) + jacobian(ExtraParts_qd,qd);%Eqn 32a    
    ExtModRNEA_qdqd_V2 = -upper;

    %Second partials of Extended Mod RNEA with q tau 
    upper =  -jacobian(ExtraParts_tau,q); %Eqn 32a 
    ExtModRNEA_qtau_V2 = upper;
  
    %Second partials of Extended Mod RNEA with q qd
    upper  = jacobian(jacobian(out_Extended+acc_act,qd),q) - jacobian(ExtraParts_qd,q); %Eqn 32a 
    ExtModRNEA_qqd_V2 = upper; 
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{SecondMat,ExtModRNEA_qtau_V2}); 
    [bc_OutSecondMat, bc_OutQTau] = second_mat_fn(q,qd,u,bc_Im_dyn,bc_Kinv_vec_T(x,bb),FOP_bc_Im_dyn,FQ_bc_Im_dyn); 
    bc_Im_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,bb},{all_bc_Im_dyn,bc_OutSecondMat,bc_OutQTau});  

%% Front Impact Dynamics

    [out_Extended,a_jnt,~,v_jnt] =modID_casadi_Extended(robot_no_grav,q,qd*0,qdd,mu,ft_force_q__MX); %H*qdd - ft_j*Fc
    out_Extended = -out_Extended + modID_casadi(robot_no_grav,q,0*qd,mu,qd); % mu'*H*qd;    
    
    acc_act = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
         
    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_q,mu,ft_q_ext_MX); 
    ftJnuq = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
    ExtraParts_q = A - ftJnuq;

    dd = jacobian(ExtraParts_q,q);
    upper = hessian(out_Extended + acc_act,q) - (dd + dd');%Eqn 32a
    
    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_qd,mu,ft_qd_ext_MX); 
    ftJnuqd = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
    ExtraParts_qd = A - ftJnuqd;  

    [A,a_jnt,~,v_jnt] = modID_casadi_Extended(robot_no_grav,q,0*qd,nu_tau,mu,ft_tau_ext_MX); 
    ftJnutau = muf'*jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
    ExtraParts_tau = A - ftJnutau;
      
    %Second partials of Extended Mod RNEA with qq
    ExtModRNEA_qq_V2 =upper;
    
    %Second partials of Extended Mod RNEA with qd
    upper = hessian(out_Extended+acc_act,qd) + jacobian(ExtraParts_qd,qd);%Eqn 32a    
    ExtModRNEA_qdqd_V2 = -upper;

    %Second partials of Extended Mod RNEA with q tau 
    upper =  -jacobian(ExtraParts_tau,q); %Eqn 32a 
    ExtModRNEA_qtau_V2 = upper;
  
    %Second partials of Extended Mod RNEA with q qd
    upper  = jacobian(jacobian(out_Extended+acc_act,qd),q) - jacobian(ExtraParts_qd,q); %Eqn 32a 
    ExtModRNEA_qqd_V2 = upper; 
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{SecondMat,ExtModRNEA_qtau_V2}); 
    [ft_OutSecondMat, ft_OutQTau] = second_mat_fn(q,qd,u,ft_Im_Dyn,ft_Kinv_vec_T(x,bb),FOP_ft_Im_dyn,FQ_ft_Im_dyn); 
    ft_Im_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,bb},{all_ft_Im_dyn,ft_OutSecondMat,ft_OutQTau}); 

%% Tensor Methods
    %% Front Stance 
      ft_xx = jacobian(jacobian(ft_dyn, x),x) ;
      ft_qu = jacobian(jacobian(ft_dyn, q),u);
      ft_tensor = Function('KKT_x', {x,u},  {all_ft_dyn,ft_xx,ft_qu} );
%       Dyn_funcs.Tens.John = Function('KKT_xx', {x,u},  {all_ft_dyn,ft_xx,jacobian(jacobian(ft_dyn, x),u)} );
    %% Back Stance 
      bc_xx = jacobian(jacobian(bc_dyn, x),x) ;
      bc_qu = jacobian(jacobian(bc_dyn, q),u);
      bc_tensor = Function('KKT_x', {x,u},  {all_bc_dyn,bc_xx,bc_qu} );

      %% free dy 
      fr_xx = jacobian(jacobian(fr_dyn, x),x) ;
      fr_qu = jacobian(jacobian(fr_dyn, q),u);
      fr_tensor = Function('KKT_x', {x,u},  {all_fr_dyn,fr_xx,fr_qu} );

      %% Front Impact Stance 
      ft_Im_xx = jacobian(jacobian(ft_Im_Dyn, x),x) ;
      ft_Im_qu = jacobian(jacobian(ft_Im_Dyn, q),u);
      ft_Im_tensor = Function('KKT_x', {x,u},  {all_ft_Im_dyn,ft_Im_xx,ft_Im_qu} );

      %% Back Impact Stance 
      bc_Im_xx = jacobian(jacobian(bc_Im_dyn, x),x) ;
      bc_Im_qu = jacobian(jacobian(bc_Im_dyn, qd),u);
      bc_Im_tensor = Function('KKT_x', {x,u},  {all_bc_Im_dyn,bc_Im_xx,bc_Im_qu} );
  
%%
% end
%{
  kkt_out = ft_dyn;
  Dyn_funcs.k_x = Function('K_x', {x,u},  {jacobian(kkt_out, x)} );
  Dyn_funcs.kkt_x = Function('KKT_x', {x,u},  { jacobian(jacobian(kkt_out, x),x) } );
  Dyn_funcs.kkt_qtau= Function('KKT_qtau', {q,qd,u},  { jacobian(jacobian(kkt_out, q),u) } );

  


  %Does it work 
  1==1;
  qVal = rand(size(q)); qdVal = rand(size(qd));
  tauVal = rand(size(u)); bbVal = rand(7,1); 
  xVal = [qVal; qdVal]; 
  
      full( Dyn_funcs.k_x(xVal,tauVal))-  full(Dyn_funcs.ft_dyn_x(xVal,tauVal))
  
  h1 = Dyn_funcs.fr_dyn_x(xVal,tauVal);
  h2 = Dyn_funcs.fr_dyn_u(xVal,tauVal);
  
  all = [h1 h2]; 
  FOPVal = all(1:Nb,:);
  fQVal = all(Nb+1:end,:);
        
  kkt_xx = full(Dyn_funcs.kkt_x(xVal,tauVal)); 
  kkt_xx = reshape(kkt_xx,[(Nb),Nb*2,Nb*2]); 
  v_xx_kkt = full(Tens3byVec(kkt_xx,bbVal','pre'))
   
  kkt_qtau = full(Dyn_funcs.kkt_qtau(qVal,qdVal,tauVal)); 
  kkt_qtau = reshape(kkt_qtau,[Nb Nb length(u)]);
  v_qtau_kkt = (full(Tens3byVec(kkt_qtau,bbVal','pre'))); 
    
  bOther = v_qtau_kkt;
  [a b] = fr_OutSecondMat_fn(xVal,tauVal,bbVal,FOPVal);
  [A B] = fr_ExtMod_OutSecondMat_fn(xVal,tauVal,bbVal,FOPVal);
    
  out = full(a - v_xx_kkt);
  out2 = full(b - bOther');
  out3 = full(A - v_xx_kkt)
  out4 = full(B - bOther');%works
    
  ad = Dyn_funcs.bc_Im_dyn(xVal);
    

    ad = rand(9,1)
    D3 = Function('dad',{q,aVec},{ft_K*aVec});
    ft_K_vec(qVal,ad) 
    D3(qVal,ad)

    D3 = Function('dad',{q,aVec},{ft_J*qdd});
    dfas = ID_casadi_Extended(robot_no_grav,q,qd*0,qdd*0,ft_force_q__MX); 
    G = Function('ad',{q,aVec},{dfas})
    G(qVal,ad)
    D3(qVal,ad)

  
 %}
    
1==1; 
%% Set All The Functions 

Dyn_funcs.Exp.Ft_stnce = ft_OutSecondMat_fn; 
Dyn_funcs.Exp.Bc_stnce = bc_OutSecondMat_fn; 
Dyn_funcs.Exp.Fr_Dyn = fr_OutSecondMat_fn; 
Dyn_funcs.Exp.Ft_Imp = ft_Im_OutSecondMat_fn; 
Dyn_funcs.Exp.Bc_Imp = bc_Im_OutSecondMat_fn;

Dyn_funcs.ExtMod.Ft_stnce = ft_ExtMod_OutSecondMat_fn; 
Dyn_funcs.ExtMod.Bc_stnce = bc_ExtMod_OutSecondMat_fn; 
Dyn_funcs.ExtMod.Fr_Dyn = fr_ExtMod_OutSecondMat_fn;
Dyn_funcs.ExtMod.Ft_Imp = ft_Im_ExtMod_OutSecondMat_fn;
Dyn_funcs.ExtMod.Bc_Imp = bc_Im_ExtMod_OutSecondMat_fn;

Dyn_funcs.Tens.Ft_stnce = ft_tensor; 
Dyn_funcs.Tens.Bc_stnce = bc_tensor; 
Dyn_funcs.Tens.Fr_Dyn = fr_tensor;
Dyn_funcs.Tens.Ft_Imp = ft_Im_tensor;
Dyn_funcs.Tens.Bc_Imp = bc_Im_tensor;

  

end
function [a_cnt,v_cnt]= jntToTask(a_jnt,v_jnt,q,robot_params,ftorbc,gravComp) 

%Front/rear foot Transformation Matrix in inertia frame
[~,T_front,T_rear] = fkin_quadruped_2D(q, robot_params); 

LinkLength = robot_params.kneelinkLength;
if ftorbc ==1 
    %what's acceleration of knee joint felt at the foot
    n = size(a_jnt{5},2);
    L = [0, 0, -LinkLength]' * ones(1,n);
    a_cnt = Vpt(a_jnt{5},L);
    v_cnt = Vpt(v_jnt{5},[0, 0, -LinkLength]');
else
   %what's acceleration of knee joint felt at the foot
   n = size(a_jnt{5},2);
   L = [0, 0, -LinkLength]' * ones(1,n); %doing this to support 3*n accelerations. Important.
   a_cnt = Vpt(a_jnt{7},L);
   v_cnt = Vpt(v_jnt{7},[0, 0, -LinkLength]');
end

%Corriolis effects of gravity : This is omega x v
if ftorbc == 1
    acc_coriolis = cross( v_jnt{5}(1:3),v_cnt);
else
    acc_coriolis = cross( v_jnt{7}(1:3),v_cnt);
end

if ftorbc ==1
    %Front foot Transformation Matrix in inertia frame
    rot = T_front(1:3,1:3); 
else
    %rear foot Transformation Matrix in inertia frame
    rot = T_rear(1:3,1:3); 
end

a_cnt = rot*(a_cnt + acc_coriolis );
v_cnt = rot*v_cnt;
if gravComp
    a_cnt = a_cnt + [0 0 -9.81]'; 
end

%Remove the y component 
a_cnt = [a_cnt(1,:);a_cnt(3,:)]; %doing this to support 3*n accelerations. Important.
v_cnt = [v_cnt(1);v_cnt(3)];

end

