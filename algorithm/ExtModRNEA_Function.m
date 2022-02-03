function Out = ExtModRNEA_Function(model,robot_params, x,u)


x_size = length(x);
u_size = length(u);
q_size = length(x)/2;
y_size = 2; % 2D quadruped


% Split the state
q = x(1:q_size,1);
qd = x(q_size+1:end,1);

%Mass matrix and Corriolis Matrix
if isa(q,'sym')
    [H, C] = HandC(model, q, qd);
else 
    import casadi.*
    [H, C] = HandC_casadi(model, q, qd);
end

%% Contact jacobian and its derivative
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

% Free Dynamics (flying-phase dynamics)
S = [zeros(q_size-u_size, u_size);
     eye(u_size)]; % selection matrix
tau_b = S*u - C; %contro, grav term, centrifugal and corriolis terms

% front stance dynamics
ft_b = [tau_b;ft_Jd*qd];
% back stance dynamics
bc_b = [tau_b;bc_Jd*qd];
%impact dynamics ; same for front and back impact dynamics
im_b = [H*qd; zeros(size(ft_Jd*qd))]; 


%Form the ft KKT matrix 
ft_K = [H, ft_J';
        -ft_J, zeros(size(ft_J,1),size(ft_J,1))]; % KKT Dynamics matrix
%Form the bc KKT matrix 
bc_K = [H, bc_J';
        -bc_J, zeros(size(bc_J,1),size(bc_J,1))]; % KKT Dynamics matrix
%%
vVec = MX.sym('vVec',[q_size+2 1]); 
ft_vVec = ft_K \ ft_b; 
bc_vVec = bc_K \ bc_b; 
fr_im_vVec = ft_K \ im_b; 
bc_im_vVec = bc_K \ im_b; 


%% First partials 

%Free Dynamics Partials
model_no_gravity = model; model_no_gravity.gravity = [0 0 0];
% Function for multiplying by M inverse (using ABA) [O(N)]
InvM_Mult = Function('InvM_Mult',{q,u}, {FDab_casadi(model_no_gravity, q, qd*0, u)}); % O(N) function
vect = MX.sym('vect',[q_size 1]); 
act_vect  = InvM_Mult(q,tau_b); % this is qdd = Hinv * b 

fr_first_q =   InvM_Mult(q,jacobian(tau_b,q)-jacobian(H*vect,q));
fr_first_qd =  InvM_Mult(q,jacobian(tau_b,qd)-jacobian(H*vect,qd));
fr_first_tau = InvM_Mult(q,jacobian(tau_b,u)); %Hinv * S 

%All Together 
fr_first_All = [fr_first_q fr_first_qd fr_first_tau]; 
%Making this function to remove 'dependency' on vect
fr_first_AllFn = Function('FR_First_All',{x,u,vect},{fr_first_All});
fr_OutFirstAll = fr_first_AllFn(x,u,act_vect); 

Out.fr_firstAll = Function('fr_firstAll',{x u},{fr_OutFirstAll});


%Front Stance Partials
% ft_first_q = Hinv_vecFn(q,qd,(jacobian(ft_b,q) - jacobian(ft_K*vVec,q)));   
% ft_first_qd = Hinv_vecFn(q,qd,jacobian(ft_b,qd));
% ft_first_tau = Hinv_vecFn(q,qd,[S;zeros(2,Nb)]); %Hinv*[S;0]

ft_first_q = ft_K \ (jacobian(ft_b,q) - jacobian(ft_K*vVec,q));   
ft_first_qd = ft_K \ (jacobian(ft_b,qd));
ft_first_tau = ft_K \ [S;zeros(2,Nb)]; %Hinv*[S;0]

%All Together 
ft_first_All = [ft_first_q ft_first_qd ft_first_tau];
%Making this function to remove 'dependency' on vVec
ft_first_AllFn = Function('FT_First_All',{x,u,vVec},{ft_first_All});
ft_OutFirstAll = ft_first_AllFn(x,u,ft_vVec); 

Out.ft_firstAll = Function('ft_firstAll',{x,u},{ft_OutFirstAll});


%Back Stance Partials
% bc_first_q = Hinv_vecFn(q,qd,(jacobian(bc_b,q) - jacobian(bc_K*vVec,q)));   
% bc_first_qd = Hinv_vecFn(q,qd,jacobian(bc_b,qd));
% bc_first_tau = Hinv_vecFn(q,qd,[S;zeros(2,Nb)]); %Hinv*[S;0]

bc_first_q = bc_K \ (jacobian(bc_b,q) - jacobian(bc_K*vVec,q));   
bc_first_qd = bc_K \ (jacobian(bc_b,qd));
bc_first_tau = bc_K \ ([S;zeros(2,Nb)]); %Hinv*[S;0]

%All Together 
bc_first_All = [bc_first_q bc_first_qd bc_first_tau];
%Making this function to remove 'dependency' on vVec
bc_first_AllFn = Function('BC_First_All',{x,u,vVec},{bc_first_All});
bc_OutFirstAll = bc_first_AllFn(x,u,bc_vVec); 

Out.bc_firstAll = Function('bc_firstAll',{x,u},{bc_OutFirstAll});


%Front Impact Dynamics 
fr_im_first_q = ft_K \ (jacobian(im_b,q) - jacobian(ft_K*vVec,q));
fr_im_first_qd = ft_K \ (jacobian(im_b,qd));
fr_im_first_tau = jacobian(im_b,u); % all zeros

%All Together 
fr_im_first_All = [fr_im_first_q fr_im_first_qd fr_im_first_tau]; 
%Making this function to remove 'dependency' on vVec 
fr_im_firstAllFn = Function('fr_im_first_all',{x u,vVec},{fr_im_first_All}); 
fr_im_OutFirstAll = fr_im_firstAllFn(x,u,fr_im_vVec); 

Out.fr_im_first_All = Function('fr_im_first_All',{x,u},{fr_im_OutFirstAll}); 


%Back Impact Dynamics 
bc_im_first_q = bc_K \ (jacobian(im_b,q) - jacobian(bc_K*vVec,q));
bc_im_first_qd = bc_K \ (jacobian(im_b,qd));
bc_im_first_tau = jacobian(im_b,u); % all zeros

%All Together 
bc_im_first_All = [bc_im_first_q bc_im_first_qd bc_im_first_tau]; 
%Making this function to remove 'dependency' on vVec 
bc_im_firstAllFn = Function('bc_im_first_all',{x u,vVec},{bc_im_first_All}); 
bc_im_OutFirstAll = bc_im_firstAllFn(x,u,bc_im_vVec); 

Out.bc_im_first_All = Function('bc_im_first_All',{x,u},{bc_im_OutFirstAll}); 


%% Second Partials




end