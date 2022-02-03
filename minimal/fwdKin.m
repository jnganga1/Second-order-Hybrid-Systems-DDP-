function [p, R, p1, R1] = fwdKin(q,model)
%SimPeriod Simulate one full period
t1 = q(1);
t2 = q(2);
t3 = q(3);
l1 = model.l1;
l2 = model.l2;


R01 = expm(cross([t1 0 0]));
R12 = expm(cross([0 t2 0]));
R23 = expm(cross([0 t3 0]));
p1 = R01*R12 * [ l1 ; 0 ;0];
R1 = R01*R12;

p = p1 + R01*R12*R23*[ l2 ; 0 ;0];
R = R01*R12*R23;


end

