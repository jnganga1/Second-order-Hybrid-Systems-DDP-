function DDP_Callback(xbar, lambda_bar, ubar, Vprev, du, K, dV, phys_params_vals,params, callback_params)

    callback_params.h_current_guess.XData = xbar(1,:);
    callback_params.h_current_guess.YData = xbar(2,:);
  
    if callback_params.plot_V_prediction

        V_store = [];
        eps_store = [];
        for eps = linspace(0,1,50)
            eps_store(end+1) = eps;
            [V, ~, ~, ~] = ForwardPass(xbar, lambda_bar, ubar,  du, K, eps, phys_params_vals,params);
            V_store(end+1) = V-Vprev;
        end
        callback_params.h_Vactual.XData = eps_store; 
        callback_params.h_Vactual.YData = V_store;

        Vpred = eps_store.*(1-eps_store/2)*dV;

        callback_params.h_Vpred.XData = eps_store; 
        callback_params.h_Vpred.YData = Vpred;

        ylim([min([V_store , Vpred] ), max([V_store, Vpred])]);
    end
    
    if callback_params.UserAdvanceIteration
        input('Continue?');
    else 
        pause(callback_params.pause);
    end
end