function [xk,g] = implicit_Nwtn(x0,u0,h,params,F,DF) 

q_size = params.q_size;
y_size = params.y_size;

xk = x0;
% xk = euler(x0,u0,h,params,F);

for i = 1:10
    KKT_eq = full(F(xk,u0)); %FD  
    qdd = KKT_eq(1:q_size,1);
    qd = xk(q_size+1:end,1); %qd_k+1
    f_con = [qd;qdd];
    G = x0 + h * f_con - xk;

    all = full(DF(xk,u0));
    qdd_x = all(1:q_size,1:q_size*2); 
    qd_x = [zeros(q_size) eye(q_size)];
    Gprime = h*[qd_x;qdd_x] - eye(q_size*2);

    if any(isnan(Gprime(:))) || cond(Gprime)> 1e6
        xk = nan*ones(q_size*2,1);
        g = nan*ones(y_size,1);
%         fprintf('Returned. h was %.f\n',h)
        return
    end
    L  = xk - Gprime\G;
%     fprintf('Error: %.4f\n',norm(L));
    if norm(xk - L ) <  1e-7    
        break 
    end 
    xk = L;
    
end


% MaxItrs = 100;
% 
% for i = 1:MaxItrs
%     KKT_eq = full(F(xk,u0)); %FD  
%     qdd = KKT_eq(1:q_size,1);
%     
%     qd = xk(q_size+1:end,1); %qd_k+1
%     
%     f_con = [qd;qdd];
%     
%     L = x0 + h*f_con;   
%     
%     if norm(xk - L) < 1e-7
% %         fprintf('Done at iter %i\n', i)
%         break 
%     end
%     xk = L;
%    
% end

g = zeros(y_size,1);
if length(KKT_eq) > q_size
    g = -KKT_eq(q_size+1:end,1);
end


end
