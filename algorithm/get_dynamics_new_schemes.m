function [f,g] = get_dynamics_new_schemes(x, u, ind, dyn, params,DynFns)
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
%             KKT_eq = full(DynFns.bc_dyn(x,u));
            F = DynFns.bc_dyn; DF = 0; 
            if strcmp(params.Itgr,params.Itgrs{2}) || strcmp(params.Itgr,params.Itgrs{4})
                DF =  DynFns.bc_dyn_derivs;
            end
            [f,g] = Intgr(x,u,params,F,DF);
        case {2,4,5}  % fl f2
%             KKT_eq = full(DynFns.fr_dyn(x,u));
            F = DynFns.fr_dyn; DF = 0; 
            if strcmp(params.Itgr,params.Itgrs{2}) || strcmp(params.Itgr,params.Itgrs{4})
                DF =  DynFns.fr_dyn_derivs;
            end           
            [f,g] = Intgr(x,u,params,F,DF);
        case 3      % fs
%             KKT_eq = full(DynFns.ft_dyn(x,u));
            F = DynFns.ft_dyn; DF = 0;
            if strcmp(params.Itgr,params.Itgrs{2}) || strcmp(params.Itgr,params.Itgrs{4})
                DF =  DynFns.ft_dyn_derivs;
            end            
            [f,g] = Intgr(x,u,params,F,DF);
    end
%     qdd = KKT_eq(1:q_size,1);
%     qd = x(q_size+1:end,1);
%     f_con = [qd;qdd];
%     f = x + f_con*dt;
%     
%     if ind == 1 || ind ==3
%         g = -KKT_eq(end-y_size+1:end,1);
%     else
%         g = zeros(y_size,1);
%     end

    
elseif strcmp(dyn,'reset')
    switch ind
        case {1,3,5}  % smooth transition
            f = x;
            g = zeros(y_size,1);
            return;
        case 2
            KKT_eq = full(DynFns.ft_Im_Dyn(x));
        case 4
            KKT_eq = full(DynFns.bc_Im_dyn(x));
    end
    %No intgr at reset map
    f = [x(1:q_size,1);
         KKT_eq(1:q_size,1)];
    g = -KKT_eq(end-y_size+1:end,1); 
end

function [xnew,gnew] = Intgr(x,u,params,F,DF)
    if strcmp(params.Itgr,params.Itgrs{1})
        %Euler Explicit Inter
        [xnew,gnew]  =  euler(x,u,params.dt,params,F);                 
     elseif strcmp(params.Itgr,params.Itgrs{2})
        %Implicit Method using Newton Method
        [xnew,gnew]  =  implicit_Nwtn(x,u,params.dt,params,F,DF);                 
    elseif strcmp(params.Itgr,params.Itgrs{3})
        %Explicit MidPoint Integration
        [xnew,gnew]  =  ExplicitMidPoint(x,u,params.dt,params,F);
    else
        %Implicit MidPoint Integration (using fixed integration)
        [xnew,gnew] = implicitMidPoint(x,u,params.dt,params,F,DF);
    end
end      
end

