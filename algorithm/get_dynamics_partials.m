function [fx,fu,gx,gu] = get_dynamics_partials(x, u, ind, dyn, params)
% ind: mode index in one gait (1,2,3,4)
% dyn: continuous dyn or reset map

%       flight      stance          impact      s2f
% f     phase dyn   phase dyn       impact dyn  I
% g     0           contact force   impulse     0
% fx    dyn par     dyn par         imp par     I
% fu    dyn par     dyn par         imp par     0
% gx    0           force par       impulse par 0
% gu    0           force par       impulse par 0

q_size = params.q_size;
x_size = params.x_size;
u_size = params.u_size;
y_size = params.y_size;
dt = params.dt;

fx = zeros(x_size, x_size);
fu = zeros(x_size, u_size);
gx = zeros(y_size, x_size);
gu = zeros(y_size, u_size);

%%
if strcmp(dyn, 'cont')
    switch ind
        case 1 % bs
            [K,b,Kx,bx,bu,qd_x] = BackStanceDyn(x, u);
            
        case {2,4} % f1 f2
            [K,b,Kx,bx,bu,qd_x] = FreeDynamics(x, u);
        case 3 % fs
            [K,b,Kx,bx,bu,qd_x] = FrontStanceDyn(x, u);    
    end
    [KKT_dyn_x, KKT_dyn_u] = get_KKT_dyn_part(K,Kx,b,bx,bu);
    
    qdd_x = KKT_dyn_x(1:q_size,:);
    f_con_x = [qd_x;qdd_x];
    f_con_u = [zeros(q_size,u_size);
        KKT_dyn_u(1:q_size,:)];
    fx = eye(x_size) + f_con_x*dt;
    fu = f_con_u*dt;
    
    if ind == 1||ind == 3 % compute GRF partial at stance mode
        gx = -KKT_dyn_x(q_size+1:end,:);
        gu = -KKT_dyn_u(q_size+1:end,:);
    end
    
elseif strcmp(dyn, 'reset')
    switch ind
        case 1 % bs to flight
            [K,b,Kx,bx,bu,~] = BackStanceDyn(x, u);
        case 2     % front impact
            [K,b,Kx,bx,bu,q_x] = FrontImpactDyn(x, u);
        case 3 % fs to flight
            [K,b,Kx,bx,bu,~] = FrontStanceDyn(x, u);
        case 4     % back impact
            [K,b,Kx,bx,bu,q_x] = BackImpactDyn(x, u);
    end
    [KKT_dyn_x, KKT_dyn_u] = get_KKT_dyn_part(K,Kx,b,bx,bu);
    
    if ind == 1|| ind ==3
        fx = eye(x_size);
    else
        fx = [q_x;
              KKT_dyn_x(1:q_size,:)];
    end
    fu = zeros(x_size, u_size); % control doesn't affect instantaneous dyn
  
    gx = -KKT_dyn_x(q_size+1:end,:);
    gu = -KKT_dyn_u(q_size+1:end,:);
end
end

%% Compute partials of KKT dynamics KKT_dyn = K\b
function [KKT_dyn_x, KKT_dyn_u] = get_KKT_dyn_part(K,Kx,b,bx,bu)
K_inv_x = zeros(size(Kx)); % partial of inv(K) w.r.t. x 
for i = 1:size(Kx,3)
    K_inv_x(:,:,i) = -K\Kx(:,:,i)/K;
end

KKT_dyn_x = zeros(size(bx)); % matrix
for i = 1:size(KKT_dyn_x,2)
    KKT_dyn_x(:,i) = K_inv_x(:,:,i)*b;    % Jacobian of KKT_dyn
end

KKT_dyn_x = KKT_dyn_x + K\bx;
KKT_dyn_u = K\bu;
end

