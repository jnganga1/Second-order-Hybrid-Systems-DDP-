function [f,g] = get_dynamics(x, u, ind, dyn, params,DynFns)
% ind: mode index in one gain
% dyn: continuous dyn or reset map
q_size = params.q_size;
y_size = params.y_size;
if size(x,1)~=params.x_size
    error('Size of x is incorrect!');
elseif size(u,1)~=params.u_size
    error('Size of u is incorrect!');
end
dt = params.dt;

if strcmp(dyn,'cont')
    switch ind
        case 1      % bs
            [K,b] = BackStanceDyn(x, u);
        case {2,4}  % fl f2
            [K,b] = FreeDynamics(x, u); % KKTx is tensor
        case 3      % fs
            [K,b] = FrontStanceDyn(x, u); 
    end
    KKT_eq = K\b; 
    
    qdd = KKT_eq(1:q_size,1);
    qd = x(q_size+1:end,1);
    f_con = [qd;qdd];
    f = x + f_con*dt;
    
    if ind == 1 || ind ==3
        g = -KKT_eq(end-y_size+1:end,1);
    else
        g = zeros(y_size,1);
    end
    
elseif strcmp(dyn,'reset')
    switch ind
        case {1,3}  % smooth transition
            f = x;
            g = zeros(y_size,1);
            return;
        case 2
            [K,b] = FrontImpactDyn(x, u);
        case 4
            [K,b] = BackImpactDyn(x, u);
    end
    KKT_eq = K\b;
    f = [x(1:q_size,1);
         KKT_eq(1:q_size,1)];
    g = -KKT_eq(end-y_size+1:end,1); 
end
        
end

