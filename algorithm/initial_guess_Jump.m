function [xbar, ybar, ubar] = initial_guess_Jump(xbar, ubar, robot_params, params,DynFns)
P = params.P;
ybar = zeros(params.y_size,params.len_sim,2);

%% flight phase reference angles
the_ref = [pi/4,-pi*7/12,pi/4,-pi*7/12]';
% the_ref = deg2rad([45,-180,45,-180]');
    
l_0 = 0.2462;
kp = 5*diag([8,1,12,10]);
kd = 1*diag([1,1,1,1]);
K_sp = 2200;
ind_bs = 0;
ind_fs = 0;

k = 1; % DDP iteration index

for mode_ind = 1:params.num_mode
    for i = 1:params.len_mode(mode_ind)
        xk = xbar(:,k);
        qk = xk(1:params.q_size,1);
        
        if i < params.len_mode(mode_ind) % continuous dyn
            switch P(mode_ind)
                case 1 % bs
                    ind_bs = ind_bs + 1;
                    if isa(xk,'sym')
                        [Jc,~] = BackJacobian(xk);
                    else
                        Jc = full(DynFns.BckJac(xk));
                    end
                    Jc(:,1:3) = [];
                    v = spring_vec(qk,robot_params, P(mode_ind));
                    Fk = -v/norm(v)*K_sp*(norm(v)-l_0);
                    uk = -Jc'*Fk*3;
                    if ind_bs ~= 1
                        uk = -Jc'*Fk*1.7;
                    end
                case {2,4} % f1 f2
                    the_k = xk(4:7,1);
                    thed_k = xk(11:end,1);
                    uk = kp*(the_ref-the_k) - kd*thed_k; % first flight
                case 3 % fs
                    ind_fs = ind_fs + 1;
                    if isa(xk,'sym')
                        [Jc,~] = FrontJacobian(xk);
                    else
                        Jc = full(DynFns.FrntJac(xk));
                    end
                    Jc(:,1:3) = [];
                    v = spring_vec(qk,robot_params, P(mode_ind));
                    Fk = -v/norm(v)*K_sp*(norm(v)-l_0);
                    uk = -Jc'*Fk*2.2;
                    if ind_fs ~= 1
                        uk = -Jc'*Fk*1.4;
                    end
            end
            [xk_next,yk] = get_dynamics_new(xk, uk, P(mode_ind),'cont',params,DynFns);
            
        else % reset map
            uk = zeros(params.u_size,1); % no control at priori instant
            [xk_next,yk] = get_dynamics_new(xk, uk, P(mode_ind),'reset',params,DynFns);
        end
        
        % store data
        ubar(:,k) = uk;
        if P(mode_ind) == 1
            ybar(:,k,2) = yk;      % bs foot GRF
        elseif P(mode_ind) == 3
            ybar(:,k,1) = yk;      % fs foot GRF
        end
        if k == params.len_sim
            break;
        end
        xbar(:,k+1) = xk_next;
        k = k + 1;
    end
end

end

