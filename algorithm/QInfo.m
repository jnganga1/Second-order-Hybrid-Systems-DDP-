function [Qx, Qu, Qxx, Quu, Qux,fgxx,fgux] = QInfo(x,x1,y, u, Vx,Vxx, ind,params,DynFns)
    [~,lx,lu,ly,lxx,luu,lux,lyy] = get_running_cost(x, y, u, ind,params);
    % 1st-order partials are all Jacobians
    
    dt = params.dt; 
    if params.iLQR
%         [fx, fu, gx, gu] = get_dynamics_partials_new(x, u, ind,'cont',params,DynFns);
        [fx, fu, gx, gu] = get_dynamics_partials_new_schemes(x,x1,u, ind,'cont',params,DynFns);
        %%
        try
            Qx = lx' + fx'*Vx + gx'*ly';
        catch 
            1==1
        end
        %%
        Qu = lu' + fu'*Vx + gu'*ly';
        Qxx = lxx + gx'*lyy*gx + fx'*Vxx*fx;
        Quu = luu + gu'*lyy*gu + fu'*Vxx*fu;
        Qux = lux + gu'*lyy*gx + fu'*Vxx*fx;
        fgxx = 0; fgux = 0;
    else
        vec = [dt*Vx(length(x)/2+1:end); -ly'];
        if params.Method == 1
            [fx, fu, gx, gu,fgxx,fgux] = get_dynamics_partials_new_2order_Exp(x, u,vec,ind,'cont',params,DynFns);
        elseif params.Method == 2
            [fx, fu, gx, gu,fgxx,fgux] = get_dynamics_partials_new_2order_ExtMod(x,u,vec,ind,'cont',params,DynFns);
        elseif params.Method == 3  
            [fx, fu, gx, gu,fgxx,fgux] = get_dynamics_partials_new_2order_Tens(x, u,vec,ind,'cont',params,DynFns);        
        end
        
        Qx = lx' + fx'*Vx + gx'*ly';
        Qu = lu' + fu'*Vx + gu'*ly';
        Qxx = lxx + gx'*lyy*gx + fx'*Vxx*fx ; %Don't forget fgxx
        Quu = luu + gu'*lyy*gu + fu'*Vxx*fu;
        Qux = lux + gu'*lyy*gx + fu'*Vxx*fx ;% Don't forget fgux
    end
    %Remove numerical errors
    Quu = 0.5*(Quu+Quu');


    
% %     if params.iLQR
%     Qx = lx' + fx'*Vx + gx'*ly';
%     Qu = lu' + fu'*Vx + gu'*ly';
%     Qxx = lxx + gx'*lyy*gx + fx'*Vxx*fx;%+ fgxx;
%     Quu = luu + gu'*lyy*gu + fu'*Vxx*fu;
%     Qux = lux + gu'*lyy*gx + fu'*Vxx*fx;%+ fgux
%     

end

