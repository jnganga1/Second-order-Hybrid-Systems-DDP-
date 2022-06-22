function [Vstore, xbar, ybar, ubar, h, hstore, Vxx,params,Vstore_no_mu] = DDP(xbar, ubar, robot_params, params,DynFns,maxItr)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initial Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vbar = 0;
du = zeros(params.u_size, params.len_sim);
K = zeros(params.u_size, params.x_size, params.len_sim);

hstore = []; 
Vstore = [];
Vstore_no_mu = []; 
[Vbar, xbar, ybar, ubar, h,V_no_mu] = forward_pass(xbar, ubar, du, K, 1, robot_params, params,DynFns);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDP Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic

iter = 1;
regularization = 0; % for full second order DDP



while 1 == 1    
    % Backward Pass 
%     if params.iLQR
% %         [dV, Vx, Vxx, du, K] = backward_pass(xbar, ybar, ubar, params,regularization,DynFns);          
%     else % It is full second order DDP, so we must be careful
        if params.Debug
            fprintf('=======================================\n');
            fprintf('Iteration %3d',iter);
        end    
        success = 0;
        while success == 0
            if params.Debug
                fprintf('\t reg=%.3e\n',regularization);
            end
            [dV, Vx, Vxx, du, K, success,regTooBig] = backward_pass(xbar, ybar, ubar, params,regularization,DynFns);
%             if params.Debug
%                 fprintf('=======================================\n');
%                 fprintf('Iteration %3d, Sucesss:%i',iter,success);
%             end            
            if success == 0
                regularization = max(regularization*4, 1e-3);
            end
            if regularization >= 10^8 || regTooBig == true
                params.DDP_OK = false;
                fprintf('Regularization exceeds boundary %e \n', 10^8);
                return;
            end
        end
%         [dV, Vx, Vxx, du, K, success] = backward_pass(xbar, ybar, ubar, params,regularization,DynFns);
        if success == 0 
            params.DDP_OK = false;
            return 
        end
        regularization = regularization / 20;
        if regularization < 1e-6
            regularization = 0;
        end
%     end
    
    % Update Plots
%     DDP_Callback(xbar, ybar, ubar, Vbar, du,K, dV ,robot_params,params, callback_params);
    Vprev = Vbar;
    
    %% Forward pass : Stepsize selection via backtracking line search
    [Vbar, xbar, ybar, ubar, h,V_no_mu] = forward_iteration(xbar, ubar, Vprev, du, K, dV , robot_params, params,DynFns);
   
    
    Change = Vprev - Vbar; 
    Vstore(end+1) = Vbar; 
    hstore(end+1) = norm(h);
    Vstore_no_mu(end+1) = V_no_mu; 
    
    iter = iter+1;
    DeltaV = params.DeltaV;
    if Change < DeltaV || iter > maxItr
        params.DDP_OK = true;
        break
    end
end
end