function [params, robot_params,Dyn_funcs] = AllTogether_set_system_functions_ALL(x, u, mode)
addpath([pwd '/support']);
addpath([pwd '/algorithm']);
import casadi.*

CodeGen = 0;

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
%% EXTENDED MODIFIED RNEA

    FOP = [nu_q nu_qd nu_tau];
    f_Q = [f_q f_qd f_tau];
    bb=[mu;muf];
%%  Front Stance Dynamics 
    mufbyX = mubyX(muf,q,robot_params,ft);
    [out_simp] = modID_casadi_Extended_ALL( robot, q, qd, qdd,mu,ft_force_q__MX,mufbyX);
    [out_qd] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_qd,mu,ft_qd_ext_MX,mufbyX);
    [out_q] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_q,mu,ft_q_ext_MX,mufbyX);              
    [out_tau] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_tau,mu,ft_tau_ext_MX,mufbyX);

    %Second partials of Extended Mod RNEA with q
    d = hessian(-out_simp,q);  b=jacobian(-out_q,q);
    ExtModRNEA_qq_V2 = d + b + b';
    
    %Second partials of Extended Mod RNEA with qd
    ExtModRNEA_qdqd_V2 = - hessian(out_simp,qd) - jacobian(out_qd,qd);
    %Second partials of Extended Mod RNEA with q tau 
    ExtModRNEA_qtau_V2 =  -jacobian(out_tau,q); 
    %Second partials of Extended Mod RNEA with q qd
    ExtModRNEA_qqd_V2  = jacobian(jacobian(-out_simp,qd),q) + jacobian(-out_qd,q); 
    
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{SecondMat,ExtModRNEA_qtau_V2}); 
    [ft_OutSecondMat, ft_OutQTau] = second_mat_fn(q,qd,u,ft_dyn,ft_Kinv_vec_T(x,bb),FOP_ft_dyn,FQ_ft_dyn); 
    ft_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,bb},{all_ft_dyn,ft_OutSecondMat,ft_OutQTau}); 
%%  Back Stance Dynamics             
    mufbyX = mubyX(muf,q,robot_params,bc);
    
    [out_simp] = modID_casadi_Extended_ALL(robot,q,qd,qdd,mu,bc_force_q__MX,mufbyX);
    
    [out_qd] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_qd,mu,bc_qd_ext_MX,mufbyX);
    [out_tau] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_tau,mu,bc_tau_ext_MX,mufbyX);
    [out_q] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_q,mu,bc_q_ext_MX,mufbyX);
    
    %kkt deriv with qd
%     Dyn_funcs.KKT_qd_only = Function('KKT_qd',{x,u},{jacobian(bc_dyn,qd)});
%     %rnea deriv with qd 
%     Dyn_funcs.Rnea_qd_only = Function('RNEA_qd',{x,u,bb},{jacobian(out_simp,qd)});%derivs of mrneac
%     %we need eeta = K^{-1}* gamma
%     HH = Dyn_funcs.Rnea_qd_only(x,u,-1.0*bc_Kinv_vec_T(x,bb)); 
%     Dyn_funcs.RNEA_qd = Function('Final',{x,u,bb},{HH});
%     %Lets evaluate and see if they're the same answer
%     
%     xVal = rand(size(x)); uVal =rand(size(u)); muVal = rand(size(mu)); mufVal = rand(size(muf)); 
%     fprintf("FO via KKT_qd");
%     full([muVal;mufVal]' * Dyn_funcs.KKT_qd_only(xVal,uVal))
%     fprintf("FO via RNEA");
%     full(Dyn_funcs.RNEA_qd(xVal,uVal,[muVal;mufVal]))
%     
%     %repeat above with q 
%     Dyn_funcs.KKT_q_only = Function('KKT_q',{x,u},{jacobian(bc_dyn,q)});
%     Dyn_funcs.Rnea_q_only = Function('RNEA_q',{x,u,aVec,bb},{jacobian(out_simp,q)}); 
%     HH = Dyn_funcs.Rnea_q_only(x,u,bc_dyn,-1.0*bc_Kinv_vec_T(x,bb)); 
%     Dyn_funcs.RNEA_q = Function('Final_q',{x,u,bb},{HH});
%     fprintf("FO via KKT_q");
%     full([muVal;mufVal]' * Dyn_funcs.KKT_q_only(xVal,uVal))
%     fprintf("FO _q via RNEA")
%     full(Dyn_funcs.RNEA_q(xVal,uVal,[muVal;mufVal]))
%     1==1; 
%     T2=jacobian(-out_q ,q);
%         Dyn_funcs.T2 = Function('T2',{x,u,FOP},{T2});

    


    %Second partials of Extended Mod RNEA with q
    d=hessian(-out_simp,q); b=jacobian(-out_q ,q);
    ExtModRNEA_qq_V2 = d + b + b';
    
    %Second partials of Extended Mod RNEA with qd
    ExtModRNEA_qdqd_V2 = -hessian(out_simp,qd) - jacobian(out_qd,qd);
    %Second partials of Extended Mod RNEA with q tau 
    ExtModRNEA_qtau_V2 =  -jacobian(out_tau,q); 
    %Second partials of Extended Mod RNEA with q qd
    ExtModRNEA_qqd_V2 = jacobian(jacobian(-out_simp,qd),q) + jacobian(-out_qd ,q); %Eqn 32a 

    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{SecondMat,ExtModRNEA_qtau_V2}); 
    [bc_OutSecondMat, bc_OutQTau] = second_mat_fn(q,qd,u,bc_dyn,bc_Kinv_vec_T(x,bb),FOP_bc_dyn,FQ_bc_dyn); 
    bc_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,bb},{all_bc_dyn,bc_OutSecondMat,bc_OutQTau}); 
%%  Free Dynamics
    out_Extended =modID_casadi(robot,q,qd,a7Vec,mu);
    ExtraParts_qd = modID_casadi(robot_no_grav,q,0*qd,mu,nu_qd)';        
    ExtraParts_q = modID_casadi(robot_no_grav,q,0*qd,mu,nu_q)';        
    ExtraParts_tau = modID_casadi(robot_no_grav,q,0*qd,mu,nu_tau)';

    %Second partials of Extended Mod RNEA with q
    d=hessian(-out_Extended,q); b=jacobian(-ExtraParts_q ,q);
    ExtModRNEA_qq_V2 = d + b + b';
    %Second partials of Extended Mod RNEA with qd
    ExtModRNEA_qdqd_V2 = -hessian(out_Extended,qd) - jacobian(ExtraParts_qd,qd);%Eqn 32a    
    %Second partials of Extended Mod RNEA with q tau 
    ExtModRNEA_qtau_V2 =  -jacobian(ExtraParts_tau,q); 
    %Second partials of Extended Mod RNEA with q qd
    ExtModRNEA_qqd_V2  = jacobian(jacobian(-out_Extended,qd),q) + jacobian(-ExtraParts_qd,q); %Eqn 32a 
    
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{x,u,a7Vec,mu,FOP},{SecondMat,ExtModRNEA_qtau_V2}); 
    [fr_OutSecondMat, fr_OutQTau] = second_mat_fn(x,u,fr_dyn,fr_Kinv_vec_T(x,mu),FOP_fr_dyn); 
    fr_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,mu},{all_fr_dyn,fr_OutSecondMat,fr_OutQTau}); 
%% Back Impact Dynamics

     mufbyX = mubyX(muf,q,robot_params,bc);
    [out_simp] = modID_casadi_Extended_ALL(robot_no_grav,q,qd*0,qdd,mu,bc_force_q__MX,mufbyX);
    out_simp = -out_simp + modID_casadi(robot_no_grav,q,0*qd,mu,qd);
    [out_qd] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_qd,mu,bc_qd_ext_MX,mufbyX);
    [out_tau] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_tau,mu,bc_tau_ext_MX,mufbyX);
    [out_q] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_q,mu,bc_q_ext_MX,mufbyX);  
    
    %Second partials of Extended Mod RNEA with qq
    b = jacobian(out_q,q);
    ExtModRNEA_qq_V2 =  hessian(out_simp,q) - (b+b');
    %Second partials of Extended Mod RNEA with qd
    ExtModRNEA_qdqd_V2= -hessian(out_simp,qd) - jacobian(out_qd,qd);%Eqn 32a    
    %Second partials of Extended Mod RNEA with q tau 
    ExtModRNEA_qtau_V2 =  -jacobian(out_tau,q); %Eqn 32a 
    %Second partials of Extended Mod RNEA with q qd
    ExtModRNEA_qqd_V2  = jacobian(jacobian(out_simp,qd),q) - jacobian(out_qd,q); %Eqn 32a 
    
    
    %All Together
    SecondMat = [ExtModRNEA_qq_V2 ExtModRNEA_qqd_V2'; 
            ExtModRNEA_qqd_V2 ExtModRNEA_qdqd_V2]; 
    second_mat_fn = Function('second_mat_fn',{q,qd,u,aVec,bb,FOP,f_Q},{SecondMat,ExtModRNEA_qtau_V2}); 
    [bc_OutSecondMat, bc_OutQTau] = second_mat_fn(q,qd,u,bc_Im_dyn,bc_Kinv_vec_T(x,bb),FOP_bc_Im_dyn,FQ_bc_Im_dyn); 
    bc_Im_ExtMod_OutSecondMat_fn = Function('ExtMod_OutSecondMat_fn',{x,u,bb},{all_bc_Im_dyn,bc_OutSecondMat,bc_OutQTau});  

%% Front Impact Dynamics
    mufbyX = mubyX(muf,q,robot_params,ft);
    [out_simp] = modID_casadi_Extended_ALL(robot_no_grav,q,qd*0,qdd,mu,ft_force_q__MX,mufbyX);
    out_simp = -out_simp + modID_casadi(robot_no_grav,q,0*qd,mu,qd);
    [out_qd] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_qd,mu,ft_qd_ext_MX,mufbyX);
    [out_tau] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_tau,mu,ft_tau_ext_MX,mufbyX);
    [out_q] = modID_casadi_Extended_ALL(robot_no_grav,q,0*qd,nu_q,mu,ft_q_ext_MX,mufbyX);  
    
    
    %Second partials of Extended Mod RNEA with qq
    b = jacobian(out_q,q);
    ExtModRNEA_qq_V2 =  hessian(out_simp,q) - (b+b');
    %Second partials of Extended Mod RNEA with qd
    ExtModRNEA_qdqd_V2= -hessian(out_simp,qd) - jacobian(out_qd,qd);%Eqn 32a    
    %Second partials of Extended Mod RNEA with q tau 
    ExtModRNEA_qtau_V2 =  -jacobian(out_tau,q); %Eqn 32a 
    %Second partials of Extended Mod RNEA with q qd
    ExtModRNEA_qqd_V2  = jacobian(jacobian(out_simp,qd),q) - jacobian(out_qd,q); %Eqn 32a 
    
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
1==1; 
%% Set All The Functions 

% Dyn_funcs.Exp.Ft_stnce = ft_OutSecondMat_fn; 
% Dyn_funcs.Exp.Bc_stnce = bc_OutSecondMat_fn; 
% Dyn_funcs.Exp.Fr_Dyn = fr_OutSecondMat_fn; 
% Dyn_funcs.Exp.Ft_Imp = ft_Im_OutSecondMat_fn; 
% Dyn_funcs.Exp.Bc_Imp = bc_Im_OutSecondMat_fn;

%Extended Mod RNEA
Dyn_funcs.ExtMod.Ft_stnce = ft_ExtMod_OutSecondMat_fn; 
Dyn_funcs.ExtMod.Bc_stnce = bc_ExtMod_OutSecondMat_fn; 
Dyn_funcs.ExtMod.Fr_Dyn = fr_ExtMod_OutSecondMat_fn;
Dyn_funcs.ExtMod.Ft_Imp = ft_Im_ExtMod_OutSecondMat_fn;
Dyn_funcs.ExtMod.Bc_Imp = bc_Im_ExtMod_OutSecondMat_fn;


%Tensor
Dyn_funcs.Tens.Ft_stnce = ft_tensor; 
Dyn_funcs.Tens.Bc_stnce = bc_tensor; 
Dyn_funcs.Tens.Fr_Dyn = fr_tensor;
Dyn_funcs.Tens.Ft_Imp = ft_Im_tensor;
Dyn_funcs.Tens.Bc_Imp = bc_Im_tensor;


if CodeGen 
    cg_options = struct;
    cg_options.casadi_real = 'real_T';
    cg_options.casadi_int = 'int_T';
    cg_options.with_header = true;
    cg_options.cpp = true; 
    
    cg = CodeGenerator('ft_ExtMod_OutSecondMat_fn',cg_options);    cg.add(ft_ExtMod_OutSecondMat_fn); cg.generate();
    cg = CodeGenerator('bc_ExtMod_OutSecondMat_fn',cg_options);    cg.add(bc_ExtMod_OutSecondMat_fn); cg.generate();
    cg = CodeGenerator('fr_ExtMod_OutSecondMat_fn',cg_options);    cg.add(fr_ExtMod_OutSecondMat_fn); cg.generate();
    cg = CodeGenerator('ft_Im_ExtMod_OutSecondMat_fn',cg_options); cg.add(ft_Im_ExtMod_OutSecondMat_fn); cg.generate();
    cg = CodeGenerator('bc_Im_ExtMod_OutSecondMat_fn',cg_options); cg.add(bc_Im_ExtMod_OutSecondMat_fn); cg.generate();
    
    cg = CodeGenerator('ft_tensor',cg_options);    cg.add(ft_tensor); cg.generate();
    cg = CodeGenerator('bc_tensor',cg_options);    cg.add(bc_tensor); cg.generate();
    cg = CodeGenerator('fr_tensor',cg_options);    cg.add(fr_tensor); cg.generate();
    cg = CodeGenerator('ft_Im_tensor',cg_options); cg.add(ft_Im_tensor); cg.generate();
    cg = CodeGenerator('bc_Im_tensor',cg_options); cg.add(bc_Im_tensor); cg.generate();
    
    cg = CodeGenerator('bc_dyn',cg_options); cg.add(Dyn_funcs.bc_dyn); cg.generate();
    cg = CodeGenerator('ft_dyn',cg_options); cg.add(Dyn_funcs.ft_dyn); cg.generate();
    cg = CodeGenerator('BckJac',cg_options); cg.add(Dyn_funcs.BckJac); cg.generate();
    cg = CodeGenerator('FrntJac',cg_options); cg.add(Dyn_funcs.FrntJac); cg.generate();
    
end


  

end

function [a_m] = mubyX(muVec,q,robot_params,ftorbc) 
%Front/rear foot Transformation Matrix in inertia frame
[~,T_front,T_rear] = fkin_quadruped_2D(q, robot_params);
a_m =  cell(length(q),1);
for i=1:length(q)
    a_m{i} = casadi.MX.zeros(1,6);
end 


B = [0,0,0,muVec(1),0,muVec(2)];
if ftorbc == 1     
  %muT * Acceleration of knee joint felt at the foot expressed in inertial Frame
  a_cnt = B * pluho(T_front);
%   a_m{5} = [a_cnt(4:end)]';
  a_m{5} = a_cnt;
else
  a_cnt = B * pluho(T_rear);
%   a_m{7} = [a_cnt(4:end)]';
   a_m{7} = a_cnt;
end

end

function [a_cnt]= CntAcc(a_jnt,v_jnt,q,robot_params,ftorbc,gravComp) 
%This Function is cleaner than jntToTask -> will clean up later 

%Front/rear foot Transformation Matrix in inertia frame
[~,T_front,T_rear] = fkin_quadruped_2D(q, robot_params); 
LinkLength = robot_params.kneelinkLength;

if ftorbc == 1     
  %Acceleration of knee joint felt at the foot expressed in inertial Frame
  a_cnt = pluho(T_front)* [0;0;0;Vpt(a_jnt{5},[0, 0, -LinkLength]')];
  corriolis = pluho(T_front)* [0;0;0;Vpt(v_jnt{5},[0, 0, -LinkLength]')];
  a_cnt = a_cnt + corriolis;
else
  a_cnt = pluho(T_rear)* [0;0;0;Vpt(a_jnt{7},[0, 0, -LinkLength]')];
  corriolis = pluho(T_rear)* [0;0;0;Vpt(v_jnt{7},[0, 0, -LinkLength]')];
  a_cnt = a_cnt + corriolis;
end

if gravComp 
    a_cnt = a_cnt - [0 0 0 0 0 9.81*1]';
end
%Remove the y component 
a_cnt = -[a_cnt(4);a_cnt(6)];

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

