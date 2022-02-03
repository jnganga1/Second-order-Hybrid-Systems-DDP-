addpath(genpath([pwd]));
addpath([pwd '\algorithm']);
addpath([pwd '\support']);
addpath([pwd '\figures']);
addpath([pwd '\output']);
addpath([pwd,'\Results']);
addpath(genpath([pwd,'\spatial_v2']));
addpath(genpath([pwd,'\spatial_v2_casadi']));
addpath('C:\Program Files (x86)\casadi-windows-matlabR2016a-v3.5.3');
import casadi.*


import casadi.*

NJ = 7; % Degree of freedom
NB = 5; % Number of free bodies
NA = 4; % Number of actuations

% state symbolic variables
q_sym = MX.sym('q',[NJ 1]); %'real'); % x, z, pitch, joint angles
q_dot_sym = MX.sym('qd', [NJ 1]);% 'real');
x_sym = [q_sym; q_dot_sym];
y_sym = MX.sym('yr',[2 1]);%'real'); % contact force (2x1) (2D quadruped)
u_sym = MX.sym('u',[NA 1]); %'real'); % actuations

load('EarthGrav_Bound')

%
params.Itgrs = {'Eulr','Imp','EulrMid','ImpMid'};
dt_lst = logspace(-7,0,10);

mass = robot_params.mass_robot; 
Inertia = robot_params.I_body(2,2); 

Mapper =  {}; 
Storage = {}; 
for iDt = 1:length(dt_lst)
    dt = dt_lst(iDt);   
    params.dt = dt; 
    
    for iItgrs = 1: length(params.Itgrs) 
        params.Itgr = params.Itgrs{iItgrs};
        [V,x,y,u,h,V_no_mu]  = forward_pass_schemes(xbar, ubar, du, K, 1, robot_params, params,DynFns);
        time = 1:length(u); time = time*params.dt;
        LinearMoment= sum(mass * x(8,:));  
        AngularMoment = sum(Inertia * x(10,:));
        Mapper{iDt,iItgrs} = ['dt:',num2str(dt),params.Itgrs{iItgrs}]; 
        Storage{iDt,iItgrs} = [LinearMoment;AngularMoment];
%         simulation(x(1:params.q_size,:),robot_params,params,false);
    end
end 

%
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

Euler = [Storage{:,1}]';
Imp  = [Storage{:,2}]';
EulrMid = [Storage{:,3}]'; 
ImpMid = [Storage{:,4}]';

iMomenta = 2;
figure; 
a1 = Euler(:,iMomenta); a1 = a1(isfinite(a1));
dt = dt_lst(isfinite(a1));
plot(dt,a1,'DisplayName','Eulr'); 
hold on
a1 = Imp(:,iMomenta); a1 = a1(isfinite(a1)); 
dt = dt_lst(isfinite(a1));
plot(dt,a1,'DisplayName','Imp') 
a1 = EulrMid(:,iMomenta); a1 = a1(isfinite(a1)); 
dt = dt_lst(isfinite(a1));
plot(dt,a1,'DisplayName','EulrMid') 
a1 = ImpMid(:,iMomenta); a1 = a1(isfinite(a1)); 
dt = dt_lst(isfinite(a1));
plot(dt,a1,'DisplayName','ImpMid') 
grid on; grid minor; legend;
L = gca; 
L.XScale ='log'
title('EARTH,angul')