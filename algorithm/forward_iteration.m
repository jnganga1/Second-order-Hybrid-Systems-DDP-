
function [V, x, y, u, h,V_no_mu] = forward_iteration(xbar, ubar, Vprev, du, K, dV, robot_params, params,DynFns)
    V = 0;
    eps = 0.05;
%     eps = 1;
    
    x = 0*xbar;
    u = 0*ubar;
    
    while eps > 1e-20
        
        % Try to take the step
        [V,x,y,u,h,V_no_mu] = forward_pass(xbar, ubar, du, K, eps, robot_params, params,DynFns);
        if params.Debug    
            fprintf('\t eps=%.3e \t DV=%.3e \t min=%.3e\n',eps, V-Vprev, params.gamma* eps*(1-eps/2)*dV );
        end
        
        % Check if V is small enough to accept the step
        if V < Vprev + params.gamma* eps*(1-eps/2)*dV
            % if so, exit
            break
        end
        
        % Else, backtrack
        eps = params.beta * eps;
    end
end
