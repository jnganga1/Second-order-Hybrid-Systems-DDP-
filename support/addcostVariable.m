function [l,l_first,l_second] = addcostVariable(l,l_first,l_second,z1,z2,delta,var,eps)
%function that takes two constraints z1 and z2 and the corresponding cost
%and its first jacobian and second jacobian and adds cost 

dt = 0.001;

% eps = 0.01;
eps = eps*dt;

%minmax control cost 
% maxU = 33;
% z1 = -var + maxu;
% z2 = var + maxu;
%minmax force normal cost 
% z1 = -var; 
% z2 = var;
%minmax cone cost
% z1 = -var; 
% z2 = var;


if length(var) == 2
  %only for the cone constraint is the jacobian not I_n
   z1_x  = [-1 1];
   z2_x  = [1 1]; 
else
    z1_x = -eye(length(var)); z1_xx = 0; 
    z2_x = eye(length(var)); z2_xx = 0;
end

%minmax barrier cost
if length(var) == 2
    [neg_cost,neg_costz,neg_costzz] =  barrier(z1,delta);
    [pos_cost,pos_costz,pos_costzz] =  barrier(z2,delta);
else 
    [neg_cost,neg_costz,neg_costzz] =  arrayfun(@barrier,z1,delta*ones(size(var)));
    [pos_cost,pos_costz,pos_costzz] =  arrayfun(@barrier,z2,delta*ones(size(var)));
end

% l = J + eb*B(z);
% lx = Jx + eb*Bz*zx  
% lxx = Jxx + eb*zx'*Bzz*zx + Bz*zxx where z_xx = 0; 
l = l+eps*sum(pos_cost+neg_cost);
l_first = l_first + eps*(pos_costz'*z2_x + neg_costz'*z1_x);
l_second = l_second + eps*(z2_x'*diag(pos_costzz)*z2_x + ...
                        z1_x' * diag(neg_costzz) * z1_x);


function [out,outz,outzz] = barrier(z,delta)
    k=2;
    if z  > delta 
        out = -log(z);
        outz = -z.^(-1);     
        outzz  = z.^(-2);
    else 
        a=(k-1)/k; 
        b=(z-k*delta)/(delta*(k-1));
        out = a *(b.^k - 1) -  log(delta);  
        outz = ((z - delta*k)/(delta*(k - 1))).^(k - 1)/delta;
        outzz = ((z - delta*k)/(delta*(k - 1))).^(k - 2);
    end
end
end