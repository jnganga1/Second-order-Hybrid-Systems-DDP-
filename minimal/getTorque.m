function [tau, q, pLeg1] = getTorque(x_c, x_d, model, leg, force)

%fprintf(1,'Leg %d\n',leg);



quat = x_c(1:4);
R1 = quatToR(quat);
pLeg = x_d.p_foot{leg} - (x_c(5:7) + R1*model.p_hip{leg});

pLeg1 = R1'*pLeg;


% hipPos = x_c(5:7)+R1*model.p_hip{leg}
% footPos= x_d.p_foot{leg}
% pLeg
% pLeg1

%pause(5);


q = invKin(pLeg1,R1,model);

t1 = q(1);
t2 = q(2);
t3 = q(3);
l1 = model.l1;
l2 = model.l2;


J = [                                      0,          - l2*sin(t2 + t3) - l1*sin(t2),         -l2*sin(t2 + t3);
     cos(t1)*(l2*sin(t2 + t3) + l1*sin(t2)),  sin(t1)*(l2*cos(t2 + t3) + l1*cos(t2)),  l2*cos(t2 + t3)*sin(t1);
     sin(t1)*(l2*sin(t2 + t3) + l1*sin(t2)), -cos(t1)*(l2*cos(t2 + t3) + l1*cos(t2)), -l2*cos(t2 + t3)*cos(t1)];
 
tau = J' * R1' * force;    

f1 = R1' * force;

tau1 = cross(pLeg1) * f1;

tau(end+1) = tau1(3);


%pause(5);
end

