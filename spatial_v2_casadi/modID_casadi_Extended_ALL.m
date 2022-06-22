function  [out] = modID_casadi_Extended_ALL( model, q, qd, qdd,lambda,f_ext,A)

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

a_grav = get_gravity(model);

import casadi.*
out = MX.zeros(size(lambda,2),1);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  wJ = S{i}*lambda(i,:);
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    w{i} = wJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i,:);
    Xa{i} = Xup{i};
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    w{i} = Xup{i}*w{model.parent(i)} + wJ;    
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i,:) + crm(v{i})*vJ;
    Xa{i} = Xup{i} * Xa{model.parent(i)};
  end
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  if length(f_ext{i}) > 0
      f{i} = f{i} - Xa{i}' \ f_ext{i};
  end
%   out = out + w{i}'*(model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i});
   out = out +  w{i}'*f{i};
% end

if nargin > 6 
    L = 0.1950;
    a_ft = cell(length(a),1);
    here = 0;
%     for i=1:length(a) 
        sz = size(a{i},2);
        a_cnt = Vpt(a{i},repmat([0 0 -L]',[1 sz]));
        a_corr = repmat(cross(v{i}(1:3),Vpt(v{i},[0 0 -L]')),[1 sz]);
        a_ft{i} = [zeros(3,sz);a_cnt + a_corr];
        here= here + A{i} * a_ft{i};
%     end
    out = out - here;
end
end
