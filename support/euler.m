function [x1,g] = euler(x0,u0,h,params,F)

q_size = params.q_size;
y_size = params.y_size;


% if strcmp(phase,'ds')
%     KKT_eq = full(DynFns.ds_dyn(x0,u0));
%     g = -KKT_eq(q_size+1:end,1);
% else %strcmp(phase,'fr')
%     KKT_eq = full(DynFns.fr_dyn(x0,u0));
%     g = zeros(y_size,1);
% end

g = zeros(y_size,1);
KKT_eq = full(F(x0,u0));
if length(KKT_eq) > q_size
    g = -KKT_eq(q_size+1:end,1);
end

qdd = KKT_eq(1:q_size,1);
qd = x0(q_size+1:end,1);
f_con = [qd;qdd];
x1 = x0 + f_con*h;


end