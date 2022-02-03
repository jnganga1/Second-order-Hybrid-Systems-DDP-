function [xk,g] = implicit_fxd(x0,u0,h,params,F) 


q_size = params.q_size;
y_size = params.y_size;

xk = x0;
% xk = euler(x0,u0,h,params,F); 
MaxItrs = 15;

for i = 1:MaxItrs
    KKT_eq = full(F(xk,u0));
    qdd = KKT_eq(1:q_size,1);
    qd = x0(q_size+1:end,1);
    f_con = [qd;qdd];
    L = x0 + h*f_con;   
    if norm(xk - L) < 1e-7
%         fprintf('Done at iter %i\n', i)
        break 
    end
    xk = L; 
end

g = zeros(y_size,1);
if length(KKT_eq) > q_size
    g = -KKT_eq(q_size+1:end,1);
end



% if strcmp(phase,'ds')
%     for i = 1:itrs
%         KKT_eq = full(DynFns.ds_dyn(xk,u0));
%         g = -KKT_eq(q_size+1:end,1);
%         qdd = KKT_eq(1:q_size,1);
%         qd = x0(q_size+1:end,1);
%         f_con = [qd;qdd];
%         l = x0 + h*f_con;
% %         fprintf('Err: %.7f\n',norm(l - xk)) 
%         xk = l; 
%     end    
% else %strcmp(phase,'fr')
%     g = zeros(y_size,1);
%     for i = 1:itrs
%         KKT_eq = full(DynFns.fr_dyn(xk,u0));
%         qdd = KKT_eq(1:q_size,1);
%         qd = x0(q_size+1:end,1);
%         f_con = [qd;qdd];
%         l = x0 + h*f_con;
% %         fprintf('Err: %.7f\n',norm(l - xk)) 
%         xk = l; 
%     end
% end





end
