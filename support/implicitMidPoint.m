function [x1,g] = implicitMidPoint(x0,u0,h,params,F,DF) 

[xhalf,~]  =  implicit_Nwtn(x0,u0,h/2,params,F,DF);
[x1,g] = euler(xhalf,u0,h/2,params,F);

end