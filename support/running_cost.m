function l = running_cost(in1,in2,in3)
%RUNNING_COST
%    L = RUNNING_COST(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    08-Feb-2022 17:19:13

q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
qd1 = in1(8,:);
qd2 = in1(9,:);
qd3 = in1(10,:);
qd4 = in1(11,:);
qd5 = in1(12,:);
qd6 = in1(13,:);
qd7 = in1(14,:);
u1 = in3(1,:);
u2 = in3(2,:);
u3 = in3(3,:);
u4 = in3(4,:);
y1 = in2(1,:);
y2 = in2(2,:);
t2 = qd1.*(2.0./2.5e1);
t3 = t2-1.0./2.5e1;
t4 = qd1-1.0./2.0;
t5 = (t3.*t4)./1.0e3;
t6 = q6./2.0e1;
t7 = q7./2.0e1;
t8 = pi.*(1.3e1./2.0e1);
t9 = q5+t8;
t10 = q3./2.0e1;
t11 = qd3.^2;
t12 = t11./1.0e5;
t13 = qd4.^2;
t14 = qd5.^2;
t15 = u1.^2;
t16 = u2.^2;
t17 = u3.^2;
t18 = u4.^2;
t19 = y2./1.0e1;
t20 = t19-1.336824;
t21 = y2-1.336824e1;
t22 = (t20.*t21)./1.0e3;
t23 = q2./2.0e2;
t24 = qd2./2.0e1;
t25 = t24-1.0./2.0e1;
t26 = qd2-1.0;
t27 = (t25.*t26)./1.0e3;
t28 = pi.*(3.0./5.0);
t29 = q5./2.0e1;
t30 = pi.*(3.0./1.0e2);
t43 = pi.*(7.0./2.0e1);
t31 = q4-t43;
t32 = q4./2.0e1;
t33 = qd6.^2;
t34 = qd7.^2;
t35 = y1.^2;
t36 = t35./1.0e4;
t37 = qd2./5.0e1;
t38 = t37+1.0./5.0e1;
t39 = qd2+1.0;
t40 = (t38.*t39)./1.0e3;
t41 = pi.*(7.0./1.0e1);
t42 = pi.*(7.0./2.0e2);
t44 = t32-pi.*(7.0./4.0e2);
t45 = (t31.*t44)./1.0e3;
t46 = pi.*(1.3e1./4.0e2);
t47 = pi./2.5e1;
t48 = q3+t47;
t49 = pi./5.0e2;
t50 = t10+t49;
t51 = (t48.*t50)./1.0e3;
t52 = t13./1.0e4;
t53 = t14./1.0e4;
t54 = t33./1.0e4;
t55 = t34./1.0e4;
t56 = t15./1.0e5;
t57 = t16./1.0e5;
t58 = t17./1.0e5;
t59 = t18./1.0e5;
l = [t5+t12+t13./5.0e3+t14./5.0e3+t15./2.0e4+t16./2.0e5+t17./1.0e6+t18./1.0e6+t22+t27+t36+t51+((q2./2.0e1+7.16e-3).*(q2+1.432e-1))./1.0e3+((q6-pi.*(7.0./2.0e1)).*(t6-pi.*(7.0./4.0e2)))./1.0e3+(t9.*(q5./5.0e1+pi.*(1.3e1./1.0e3)))./1.0e3+(t31.*(q4./5.0e1-pi.*(7.0./1.0e3)))./1.0e3+((q7+t28).*(t7+t30))./1.0e3;t5+t12+t40+t52+t53+t54+t55+t56+t57+t58+t59+((q2+1.791e-1).*(t23+8.955e-4))./1.0e3+(t9.*(t29+t46))./1.0e3+((q6-pi.*(3.0./1.0e1)).*(t6-pi.*(3.0./2.0e2)))./1.0e3+((q3-pi./3.5e1).*(t10-pi./7.0e2))./1.0e3+((q4-pi.*(1.1e1./5.0e1)).*(t32-pi.*(1.1e1./1.0e3)))./1.0e3+((q7+t41).*(t7+t42))./1.0e3;t5+t12+t15./1.0e6+t16./1.0e6+t17./2.0e4+t18./2.0e4+t22+t27+t33./5.0e3+t34./5.0e3+t36+t45+((q2+1.831e-1).*(t23+9.155e-4))./1.0e3+((q3+pi./4.0e1).*(t10+pi./8.0e2))./1.0e3+((q7+pi.*(3.0./4.0)).*(q7./5.0e1+pi.*(3.0./2.0e2)))./1.0e3+((q6-pi.*(3.3e1./1.0e2)).*(q6./5.0e1-pi.*6.6e-3))./1.0e3+((q5+t28).*(t29+t30))./1.0e3;t5+t12+t40+t45+t51+t52+t53+t54+t55+t56+t57+t58+t59+((q2+1.785e-1).*(t23+8.925e-4))./1.0e3+((q6-pi./4.0).*(t6-pi./8.0e1))./1.0e3+((q7+t8).*(t7+t46))./1.0e3+((q5+t41).*(t29+t42))./1.0e3];
