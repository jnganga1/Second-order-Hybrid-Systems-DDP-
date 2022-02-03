function [l,lx,lu,ly,lxx,luu,lux,lyy] = get_running_cost(x,y,u,ind,params)
% ind indicates which mode (in one gait cycle) is active
delta = params.delta; 
coeff = params.coeff;
eps = params.eps;
epsU = params.epsU;
l_mat = running_cost(x, y, u);
l = l_mat(ind);
1==1;

if nargout > 1
    [lx_mat,lu_mat,luu_mat,lxx_mat,lux_mat,ly_mat,lyy_mat] = lInfo(x, y, u);
    lx = lx_mat(ind,:);
    lu = lu_mat(ind,:);
    ly = ly_mat(ind,:);
    
    luu = squeeze(luu_mat(:,:,ind));
    lxx = squeeze(lxx_mat(:,:,ind));
    lux = squeeze(lux_mat(:,:,ind));
    lyy = squeeze(lyy_mat(:,:,ind));
    
    %Control cost
    maxu = 50; 
    z1 = -u + maxu; z2 = u + maxu; 
    [l,lu,luu] = addcostVariable(l,lu,luu,z1,z2,delta,u,epsU);
%     [l,lu,luu] = addcost(l,lu,luu,delta,u); 

    %minmax force normal cost 
    z1 =  - y(2); z2 = y(2); 
    [l,ly,lyy] = addcostVariable(l,ly,lyy,z1,z2,delta,y(2),eps);
    %friction cone constraint cost 
    z1 = -y(1) + coeff* y(2); z2 =  y(1) + coeff* y(2); 
    [l,ly,lyy] = addcostVariable(l,ly,lyy,z1,z2,delta,y,eps);
%     %

    

end
end

