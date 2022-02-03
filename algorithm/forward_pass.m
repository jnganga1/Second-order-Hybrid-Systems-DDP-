function [V, x, y, u, h,V_no_mu] = forward_pass(xbar, ubar, du, K, eps, robot_params, params,DynFns)
    P = params.P;
    
    V = 0;
    x = 0*xbar;
    u = 0*ubar;
    y = zeros(params.y_size, params.len_sim,2);        % 1 front 2 back
    h = zeros(1,params.num_tdconstr);                     % constraint violation
    mu = [0,0];                                        % intialize AL params
    
    x(:,1) = xbar(:,1);
    xk = xbar(:,1);
    k = 1;
    ind_tdconstr = 0;
    delta = params.delta;
    
    for mode_ind = 1:params.num_mode
        for i = 1:params.len_mode(mode_ind)
            xk = x(:,k);
            dxk = xk - xbar(:,k);
            % continuous dyn
            if i < params.len_mode(mode_ind)
                % Update with stepsize and feedback
                uk = ubar(:,k) + eps*du(:,k) + K(:,:,k)*dxk;
%                 [xk_next,yk] = get_dynamics(xk, uk, P(mode_ind),'cont',params);
                [xk_next,yk] = get_dynamics_new_schemes(xk, uk, P(mode_ind),'cont',params,DynFns);
                
                % update cost
                V = V + get_running_cost(xk,yk,uk,P(mode_ind),params);
            else
                uk = zeros(params.u_size,1);
%                 [xk_next,yk] = get_dynamics(xk, uk, P(mode_ind),'reset',params);
                [xk_next,yk] = get_dynamics_new_schemes(xk, uk, P(mode_ind),'reset',params,DynFns);

                % update cost
                if P(mode_ind)==2 || P(mode_ind)==4 || P(mode_ind)==5
                    ind_tdconstr = ind_tdconstr + 1;
                    mu = [params.mu1(ind_tdconstr),params.mu2];    % update AL params
                    h(ind_tdconstr) = violation_meas(xk, robot_params, P(mode_ind));
                end
                V_no_mu = V +  get_final_cost_no_mu(xk,P(mode_ind));
                V = V + get_final_cost(xk,mu,P(mode_ind));
            end
            
            % store data
            if P(mode_ind) == 1
                y(:,k,2) = yk;      % bs foot GRF
            elseif P(mode_ind) == 3
                y(:,k,1) = yk;      % fs foot GRF
            end
            if k == params.len_sim
                break;
            end
            u(:,k) = uk;
            try
                x(:,k+1) = xk_next;
            catch 
                1==1; 
            end
            k = k + 1;
        end
    end
    
    
end