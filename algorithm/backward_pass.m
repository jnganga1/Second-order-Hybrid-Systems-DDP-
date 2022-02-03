function [dV, Vx, Vxx, du, K, success,regTooBig] = backward_pass(xbar, ybar, ubar, params, regularization,DynFns)
    
    success = 1;
    P = params.P;
    
    largest = 0; 
    regTooBig = false;


    % Initialization
    Vx = zeros(params.x_size, params.len_sim);
    Vxx = zeros(params.x_size, params.x_size, params.len_sim);
    
    du = zeros(params.u_size, params.len_sim);
    
    K = zeros(params.u_size, params.x_size, params.len_sim);
        
    xf = xbar(:,params.len_sim);
    
    mu = [params.mu1(end),params.mu2];                     % initialize AL params
    [~,Vx(:   ,params.len_sim), Vxx(:,:,params.len_sim)] = get_final_cost(xf, mu, P(end)); 
    dV = 0;   
    dt = params.dt;
    
    ind_tdconstr = params.num_tdconstr;   % count the impact events backward
    k = params.len_sim - 1;         % back propagate the time
    
    for mode_ind = params.num_mode:-1:1
        if mode_ind == params.num_mode
            len_mode_act = params.len_mode(end)-1;      % length of current active mode
        else
            len_mode_act = params.len_mode(mode_ind);
        end
        for i = len_mode_act:-1:1
            if strcmp(params.Itgr,params.Itgrs{1})
                %Euler Explicit
                xk  = xbar(:,k);
                xk1 = xk; 
            else 
                %Implicit Euler and all other Itgrs
                xk = xbar(:,k);
                xk1 = xbar(:,k+1);
            end
            uk = ubar(:,k);
            Vxk = Vx(:,k+1);
            Vxxk = Vxx(:,:,k+1);
            
            
            if i == params.len_mode(mode_ind) % reset maps 
                if P(mode_ind)== 2 || P(mode_ind)==4 || P(mode_ind)==5
%                     ind_tdconstr = ind_tdconstr - 1;
                    mu = [params.mu1(ind_tdconstr),params.mu2];
                end
                
                [~,Phi_x,Phi_xx] = get_final_cost(xk, mu, P(mode_ind)); % jacobian and hessian of terminal cost
                if params.iLQR                                  
                    [Pr_x, ~, ~, ~] = get_dynamics_partials_new_schemes(xk,xk1,uk, P(mode_ind),'reset',params,DynFns); % jacobian of reset map
                    Vx(:,k)    =  Phi_x + (Vxk'*Pr_x)';
                    Vxx(:,:,k) =  Phi_xx + (Pr_x'*Vxxk*Pr_x)';
                else
                    vec = [Vxk(length(xk)/2+1:end)*dt; 0;0];
                    if params.Method == 1
                        [Pr_x, ~, ~, ~,fgxx,~] = get_dynamics_partials_new_2order_Exp(xk,uk,vec, P(mode_ind),'reset',params,DynFns);
                    elseif params.Method == 2
                        [Pr_x, ~, ~, ~,fgxx,~] = get_dynamics_partials_new_2order_ExtMod(xk,uk,vec,P(mode_ind),'reset',params,DynFns);
                    elseif params.Method == 3  
                        [Pr_x, ~, ~, ~,fgxx,~] = get_dynamics_partials_new_2order_Tens(xk,uk,vec,P(mode_ind),'reset',params,DynFns);        
                    end                 
                    Vx(:,k)    =  Phi_x + (Vxk'*Pr_x)';
                    Vxx(:,:,k) =   Phi_xx + (Pr_x'*Vxxk*Pr_x)'+ fgxx;
                end

                
%                 Vx(:,k)    =  Phi_x + (Vxk'*Pr_x)';
%                 Vxx(:,:,k) =  Phi_xx + (Pr_x'*Vxxk*Pr_x)';
                du(:,k) = zeros(params.u_size,1); % no control update across the transition point
                K(:,:,k) = K(:,:,k+1);
                
            elseif i < params.len_mode(mode_ind) % running
                if P(mode_ind)==1 % bs
                    yk = ybar(:,k,2);
                elseif P(mode_ind)==3 % fs
                    yk = ybar(:,k,1);
                else % flight
                    yk = zeros(params.y_size,1);
                end
                [Qx, Qu, Qxx, Quu, Qux,fgxx,fgux] = QInfo(xk,xk1,yk, uk, Vxk, Vxxk, P(mode_ind),params,DynFns);
                if params.iLQR == 0
                    % If we are doing full second order DDP, add in regularization  
                    if params.regularizationMethod == 1 
                            Qxx = Qxx + fgxx; 
                            Qux = Qux + fgux; 
                            Qxx = 0.5 *(Qxx + Qxx');
                            %{                    %Traditional Way to do it 
                            Qxx = Qxx + eye(params.x_size)*regularization;                    
                            Quu = Quu + eye(params.u_size)*regularization; 
                            % Make sure Quu is PD, if not, exit and increase regularization
%                             H = [Qxx+1e4*eye(length(Qxx)) Qux'; Qux Quu];
%                             r = length(H);
%                             [~, p] = chol(H-eye(r)*1e-9);
                            [~, p] = chol(Quu-eye(params.u_size)*1e-9);
                            if p ~= 0
                                success = 0;
                                return
                            end
                    elseif params.regularizationMethod == 4 
                            Qxx = Qxx + fgxx; 
                            Qux = Qux + fgux; 
                            Qxx = 0.5 *(Qxx + Qxx');                            
                            % Not Quite Traditional Way to do it 
                            Qxx = Qxx + eye(params.x_size)*regularization;                    
                            Quu = Quu + eye(params.u_size)*regularization; 
                            % Make sure Quu is PD, if not, exit and increase regularization
                            H = [Qxx Qux'; Qux Quu];
                            r = length(H);
                            [~, p] = chol(H-eye(r)*1e-9);
                            if p ~= 0
                                success = 0;
                                return
                            end
                    else
                        [Qxx,Qux,Quu] =regularizer_module(Qxx,Qux,Quu,fgxx,fgux,params);                        
                    end
                end
                
                % Standard equations
                Quu_inv = eye(params.u_size)/Quu;
                Quu_inv = (Quu_inv + Quu_inv')/2;
                du(:,k)  = -Quu_inv*Qu;
                K(:,:,k) = -Quu_inv*Qux;
                Vx(:,k)    = Qx  - (Qux')*(Quu_inv* Qu );
                Vxx(:,:,k) = Qxx - (Qux')*(Quu_inv* Qux);
                dV = dV - Qu'*(Quu\Qu);
            end
            if k == 1
                break
            end
            k = k - 1;
        end
    end
    
   
end