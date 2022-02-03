function [x1,g] = ExplicitMidPoint(x0,u0,h,params,F)

[xhalf,~]  =  euler(x0,u0,h/2,params,F);
[x1,g] = euler(xhalf,u0,h/2,params,F);
   

end