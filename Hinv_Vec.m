function [ft_out] = Hinv_Vec(x,model,robot_params,ft_ext,bc_ext)

import casadi.*
Nb = model.NB;    
q_size = length(x)/2; Nb=q_size;
% Casadi symbolics
% Split the state
q = x(1:q_size,1);
qd = x(q_size+1:end,1);
qdd = MX.sym('qdd',[Nb 1]);
tau = MX.sym('tau',[Nb 1]);
y1 =  MX.sym('y1',[Nb 1]);
y2 = MX.sym('y2',[2 1]);
Fc_MX =MX.sym('Fc_MX',[2 1]);
%%
%Model with zero grav
model_no_gravity = model;
model_no_gravity.gravity = [0;0;0]';

% Function for multiplying by M inverse (using ABA) [O(N)]
InvM_Mult = Function('InvM_Mult',{q,tau}, {FDab_casadi(model_no_gravity, q, qd*0, tau)}); % O(N) function

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

%% Step 1
% Front 

[~,a_jnt,~, v_jnt]= ID_casadi_Extended(model,q,qd*0,qdd,ft_ext);        

[ft_a_out,ft_v_out]=jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
as=Function('ad',{q,qd,qdd},{ft_a_out});
as2=Function('ad2',{q,qd,qdd},{ft_J*qdd + ft_Jd*qd});

as3=Function('ad3',{q,qd,qdd},{ft_v_out});
as4=Function('ad4',{q,qd,qdd},{ft_J*qd});

qVal = rand(size(q)); qdVal = rand(size(qd)); qddVal = rand(size(qdd)); 
as(qVal,qdVal,qddVal)
as2(qVal,qdVal,qddVal)
as3(qVal,qdVal,qddVal)
as4(qVal,qdVal,qddVal)

% Back 
[~,a_jnt,Xup, v_jnt]= ID_casadi_Extended(model,q,qd,qdd,bc_ext);        

[bc_a_out,bc_v_out]=jntToTask(a_jnt,v_jnt,q,robot_params,bc,1)
as=Function('ad',{q,qd,qdd},{bc_a_out});
as2=Function('ad2',{q,qd,qdd},{bc_J*qdd + bc_Jd*qd});

as3=Function('ad3',{q,qd,qdd},{bc_v_out});
as4=Function('ad4',{q,qd,qdd},{bc_J*qd});

qVal = rand(size(q)); qdVal = rand(size(qd)); qddVal = rand(size(qdd)); 
as(qVal,qdVal,qddVal)
as2(qVal,qdVal,qddVal)
as3(qVal,qdVal,qddVal)
as4(qVal,qdVal,qddVal)

% ign = cell([1 Nb]);
% [~,a_jnt,~, v_jnt]= ID_casadi_Extended(model,q,qd*0,qdd,ign);        
% [ft_a_out,ft_v_out]=jntToTask(a_jnt,v_jnt,q,robot_params,ft,0);
% 
% 
% 

1==1;

%% step 2 

% m = 2; 
% Int_inv = MX.zeros(m,m);
% 
% ei_cell = cell([1 Nb]); 
% for j=1:m
%     ei = MX.zeros(m,1); ei(j,1) = -1;
%     ei_cell{Nb} = Fpt(ei,[x_pos;y_pos]); %Format required by spatialV2
% %     ei_cell{Nb} = [0;ei];
%     [~, a ,v, ~] = FDab_casadi(model_no_gravity, q, qd*0,tau*0,ei_cell);
%      acc_act= Vpt(a{end},[1/Nb;0]);
%      vel_act = Vpt(v{end},[1/Nb;0]);
%      a = cross( [0 0 v{end}(1)]',[vel_act ; 0]);
%      acc_coriolis = cross( [0 0 v{end}(1)]',[vel_act ; 0]);
%      acc_act = R_end'*(acc_act+acc_coriolis(1:2)); 
%      1==1;
%      Int_inv(:,j)= acc_act; 
%     
% end



% qVal = rand(Nb,1); 
% Int_Inv_Val = Function('Int_Inv_Val',{q},{Int_inv});
% Int_Inv_Val(qVal)
% JJ = Function('InertiaC_inv',{q},{InertiaC_inv});
% JJ(qVal)


%% Other

m = 2; 
ft_Int_inv = MX.zeros(m,m);
ft_ei_cell = cell([1 Nb]);

bc_Int_inv = MX.zeros(m,m);
bc_ei_cell = cell([1 Nb]);

for j=1:m
    ei = MX.zeros(3,1); ei(j,1) = -1;
    ft_ei_cell{5} = Fpt(ei,Pos_ctact_3d(:,ft)); %Format required by spatialV2
    [ft_qdd, a ,v, ~] = FDab_casadi(model_no_gravity, q, qd*0,tau*0,ft_ei_cell);
    [a_cnt,~] = jntToTask(a,v,q,robot_params,ft,0);
    ft_Int_inv(:,j)= a_cnt; 
    
    bc_ei_cell{7} = Fpt(ei,Pos_ctact_3d(:,bc)); %Format required by spatialV2
    [bc_qdd, a ,v, ~] = FDab_casadi(model_no_gravity, q, qd*0,tau*0,bc_ei_cell);
    [a_cnt,~] = jntToTask(a,v,q,robot_params,bc,0);
    bc_Int_inv(:,j)= a_cnt; 

end

%The formula for contact inertia and its inverse
InertiaC = ft_J*InvM_Mult(q,ft_J'); 
ft_InertiaC_inv = InertiaC \ MX.eye(length(InertiaC)); %Always 2*2. Easy.

%The formula for contact inertia and its inverse
bc_InertiaC = bc_J*InvM_Mult(q,bc_J'); 
bc_InertiaC_inv = bc_InertiaC \ MX.eye(length(bc_InertiaC)); %Always 2*2. Easy.

%Front
fprintf('Compare contact inertia for front\n')
ft_Int_Inv_Val = Function('Int_Inv_Val',{q,qd,qdd},{ft_Int_inv});
ft_Int_Inv_Val = Function('Int_Inv_Val', {q,qd},{ft_Int_Inv_Val(q,qd,ft_qdd)});
ft_Int_Inv_Val(qVal,qdVal) 

JJ = Function('InertiaC_inv',{q},{ft_InertiaC_inv});
JJ(qVal)

%Back 
fprintf('Compare contact inertia for back\n')
bc_Int_Inv_Val = Function('Int_Inv_Val',{q,qd,qdd},{bc_Int_inv});
bc_Int_Inv_Val = Function('Int_Inv_Val', {q,qd},{bc_Int_Inv_Val(q,qd,bc_qdd)});
bc_Int_Inv_Val(qVal,qdVal) 

JJ = Function('InertiaC_inv',{q},{bc_InertiaC_inv});
JJ(qVal)

%Till you fix the loop, just use the correct stuff
ft_Int_inv = ft_InertiaC_inv;
bc_Int_inv = bc_InertiaC_inv;




%% Step 3  This is Good
% y1 = MX.sym('y1',[Nb 1]);
a1 = ft_J*InvM_Mult(q,y1);
y = [y1;y2];

x2_v1 = -ft_Int_inv * (a1 - y2); % Dim of Int_inv is wrong, I think 
% x2_v2 = -ft_InertiaC_inv * (a1-y2);
% x2_v1F = Function('x2_v1',{q,y},{x2_v1})
% x2_v2F = Function('x2_v2',{q,y},{x2_v2})
% 
% qVal = rand(Nb,1); qdVal = rand(Nb,1); tauVal = rand(Nb,1); 
% y1Val = rand(Nb,1); y2Val = rand(2,1)
% yVal = [y1Val;y2Val];
% fprintf('\nDoes x2 work\n')
% x2_v1F(qVal,yVal)
% x2_v2F(qVal,yVal)


%% Step  4 This is Good
F_fake = cell([1 Nb]); 
x2_v1_fk = [x2_v1(1,:); 0*x2_v1(1,:) ; x2_v1(2,:)]; 
F_fake{5} = Fpt(x2_v1_fk,Pos_ctact_3d(:,ft));

x1 = FDab_casadi(model_no_gravity, q, qd*0, y1,F_fake);
% x1_f = Function('x1',{q,y1,y2}, {x1});

ft_out = Function('Out',{x,y},{[x1;x2_v1]});  

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

