
function [Storage] = DDP_2DQuadruped(iLQR,Method,regularizationMethod,Itgr,OtherTests)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iLQR == 0 refers to DDP and 1 refers to iLQR 
% Method: Applies to DDP alone
% Method == 1 refers Explicit Deriv Calculation, 2 is Extended Modified RNEA, 3 is Tensor
% regularizationMethod: Applies to DDP alone 
% regularizationMethod == 1 refers to traditional regularization, i.e., repeats backward pass
% regularizationMethod == 4 refers to traditional regularization for H=[Qxx Qux';Qux Quu]
% regularizationMethod == 2 refers to point-by-point regularization 
% regularizationMethod == 3 refers to pullback of second-order terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%
% iLQR = 0; 
% Method = 1; %1 is Exp, 2 is ExtMod, 3 is Tensor

NJ = 7; % Degree of freedom
NB = 5; % Number of free bodies
NA = 4; % Number of actuations

% state symbolic variables
q_sym = MX.sym('q',[NJ 1]); %'real'); % x, z, pitch, joint angles
q_dot_sym = MX.sym('qd', [NJ 1]);% 'real');
x_sym = [q_sym; q_dot_sym];
y_sym = MX.sym('yr',[2 1]);%'real'); % contact force (2x1) (2D quadruped)
u_sym = MX.sym('u',[NA 1]); %'real'); % actuations



% set up dynamics and retrieve robot_params
% [~, ~,DynFns] = AllTogether_set_system_functions_ALL(x_sym, u_sym, '2D');
[params, robot_params,DynFns] = AllTogether_set_system_functions_ALL(x_sym, u_sym, '2D');


% state initialization

vd=0.5; 
% if OtherTests.randomInit == true
%     vd = vd + OtherTests.rndValue;
% end
q0 = [0,-0.1485,-pi/25,0.35*pi,-0.7*pi,0.25*pi,-0.65*pi]';
% q0 = [0,-0.1085,-pi/75,0.35*pi,-0.7*pi,0.35*pi,-0.65*pi]';
qd0 = [vd,-0.2,0,0,0,0,0]';
x0 = [q0;qd0];
u0 = zeros(NA,1);

% PosturePlot(q0, robot_params)

[params, callback_params] = set_sim_params(params);
params.G = [0 -robot_params.mass_robot*robot_params.g]';  
params.DeltaV = 1e-9;

%% Set which method to use and if iLQR

params.Method = Method; 
params.iLQR = iLQR;
params.regularizationMethod = regularizationMethod; 


params.Itgrs = {'Eulr','Imp','EulrMid','ImpMid'};
params.Itgr = params.Itgrs{Itgr};

%%

% Cost to go setup :: 
set_cost(x_sym, u_sym, y_sym, robot_params, params);

% Initialize
ubar = zeros(params.u_size, params.len_sim);
xbar = zeros(params.x_size, params.len_sim);
ybar = zeros(params.y_size, params.len_sim, 2); % 1 front 2 back
xbar(:,1) = x0;
Vstore =[]; Vstore_no_mu=[]; Hstore = [];
% initial guess
if OtherTests.Run ~= true
    [xbar, ybar, ubar] = initial_guess(xbar, ubar, robot_params, params,DynFns);
    simulation(xbar(1:params.q_size,:),robot_params,params,false);
end
%%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Augmented Lagrangian (AL) Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.mu1 = zeros(params.num_tdconstr,1); % Lagrangian Multipliers
params.mu2 = 5;                         % penalty multiplier
Beta_AL = 8;                            % penalty update coefficient
params.delta = 10;
params.coeff = 3;
params.umax = 34;
% params.eps=0; 
% params.epsU = 0;
params.eps =4; 
params.epsU = 100;

hbar = zeros(params.num_tdconstr,5);      % initialize violation measure for 20 AL iterations
norm_hbar = zeros(1,20);
% an = animatedline;
k_AL = 1;                               % outer loop (AL) index

du = zeros(params.u_size, params.len_sim);
K = zeros(params.u_size, params.x_size, params.len_sim);
[~, ~, ~, ~, h0,~] = forward_pass(xbar, ubar, du, K, 1, robot_params, params,DynFns);
h0 = norm(h0);
%
%% DDP for some iterations to get certain aspects 
if OtherTests.Run == true
    params.dt = OtherTests.dt; 
    maxItr = OtherTests.maxItr; 
    
    %Vary control with noise 
    if OtherTests.varyU 
%        data= load('crcMatData_New/Trad_Methods_compare_WithReg.mat');
       if OtherTests.BeforeNoise 
            data = load('MatData/VaryMethods_Trad_Data.mat');
            if OtherTests.iLQR
                xbar = data.iLQR_store.xbar{1}; 
                ubar = data.iLQR_store.ubar{1}; 
                robot_params = data.iLQR_store.rbtparams; 
                params = data.iLQR_store.params;
            else
                xbar = data.DDP_ExtMod_store.xbar{1}; 
                ubar = data.DDP_ExtMod_store.ubar{1}; 
                robot_params = data.DDP_ExtMod_store.rbtparams; 
                params = data.DDP_ExtMod_store.params;
            end
       else
            data = load('MatData/FullRuns_All_Methods_NoisedCntrlStates_BeforeNoise.mat');
            xbar = data.ExtMod_Nsd_nowOpt.xbar{1}; 
            ubar = data.ExtMod_Nsd_nowOpt.ubar{1}; 
            robot_params = data.ExtMod_Nsd_nowOpt.rbtparams; 
            params = data.ExtMod_Nsd_nowOpt.params;
        end
            
        params.Debug = 1;
        params.DeltaV = 1e-9;
        
        params.iLQR = OtherTests.iLQR; 
        params.regularizationMethod = OtherTests.RegMethod;
        params.Method = OtherTests.Method;
        
        noise = OtherTests.noise;
        StateNoise = OtherTests.StateNoise;
        
        dt = params.dt;
        t = 0:dt:(length(ubar)-1)*dt;
        figure; 
        g={};
        for i = 1:4
            subplot(2,2,i); hold on;
            g{1}=plot(t,ubar(i,:),'LineWidth',2,'DisplayName','u'); 
        end
        ubar = ubar + noise(:,1:length(ubar));
        for i = 1:4
            subplot(2,2,i); hold on;
            g{2}=plot(t,ubar(i,:),'LineWidth',2,'DisplayName','u + noise'); 
            grid on; grid minor;
            h = get(gca,'Children');
            set(gca,'Children',[flipud(h)]); 
            ylabel(['u ',num2str(i)]);
            xlabel('Time','Interpreter','latex');
        end
        legend([g{:}]);
        sgtitle(['Multiplicative = ' num2str(OtherTests.scalarNse)],'Interpreter','latex')
        savefig('Analysis/ControlNoise.fig');
        %Plot States
        figure; G = {};
        for i = 2:7
            subplot(2,3,i-1); hold on;
            G{1}=plot(t,xbar(i,:),'LineWidth',2,'DisplayName','x'); 
            grid on; grid minor;
            ylabel(['x',num2str(i)],'Interpreter','latex');
            xlabel('Time','Interpreter','latex');
        end
        x0 = xbar(:,1); 
        xbar = xbar + StateNoise(:,1:length(xbar)); 
        xbar(:,1) =x0;
        for i = 2:7
            subplot(2,3,i-1); hold on;
            G{2}=plot(t,xbar(i,:),'LineWidth',2,'DisplayName','x+Noise'); 
            grid on; grid minor;
            h = get(gca,'Children');
            ylabel(['x',num2str(i)],'Interpreter','latex');
            xlabel('Time','Interpreter','latex');
            set(gca,'Children',[flipud(h)]);            
        end
        legend([G{:}]);
        sgtitle(['Multiplicative = ' num2str(OtherTests.StateScalar)],'Interpreter','latex')
        savefig('Analysis/StateNoise.fig');
        params.filename = 'Analysis/StateNoise.gif';
        simulation(xbar(1:params.q_size,:),robot_params,params,false);

        1==1; close all;

    else 
        %VaryV Initial Condition
        data= load('MatData/Trad_Methods_compare.mat');
        xbar = data.iLQR_store.xbar{5}; 
        ubar = data.iLQR_store.ubar{5}; 
        robot_params = data.iLQR_store.rbtparams; 
        params = data.iLQR_store.params;
        params.Debug = 1;
        
        params.iLQR = OtherTests.iLQR; 
        params.regularizationMethod = OtherTests.RegMethod;
        params.Method = OtherTests.Method;
        
        xbar(8,1) = OtherTests.vd;
    end    
    
    Storage.OrigOptim = data; 
    [V0, ~, ~, ~, h0,~] = forward_pass(xbar, ubar, du, K, 1, robot_params, params,DynFns);
    start = tic;  
    [Vbar, xbar, ybar, ubar, h, hstore, Vxx,params,Vbar_no_mu] = DDP(xbar, ubar, robot_params, params,DynFns,maxItr);
    Storage.Time = toc(start);
    Storage.Vbar{k_AL} = Vbar; 
    Storage.VnoMu{k_AL} = Vbar_no_mu;
    Storage.xbar{k_AL} = xbar; 
    Storage.ybar{k_AL} = ybar; 
    Storage.ubar{k_AL} = ubar; 
    Storage.h{k_AL} = h; 
    Storage.hstore{k_AL} = hstore; 
    Storage.V0 = V0;
    Storage.OtherTests = OtherTests;
    
    Storage.rbtparams = robot_params; 
    Storage.params = params;  
    Storage.success = params.DDP_OK;  
    Storage.h0 = h0;
    fprintf('\nHouston, you have control now.\n');
    return
end
%% Full DDP code 

maxItr = 100; 
start =tic;
while 1==1 
    fprintf('=======================================\n');
    fprintf('Outer loop iteration %3d\n',k_AL);
%     if k_AL == 1 
%         maxItr = 200;
%     end

    [Vbar, xbar, ybar, ubar, h, hstore, Vxx,params,Vbar_no_mu] = DDP(xbar, ubar, robot_params, params,DynFns,maxItr);
    Storage.Vbar{k_AL} = Vbar; 
    Storage.VnoMu{k_AL} = Vbar_no_mu;
    Storage.xbar{k_AL} = xbar; 
    Storage.ybar{k_AL} = ybar; 
    Storage.h{k_AL} = h; 
    Storage.ubar{k_AL} = ubar; 
    Storage.hstore{k_AL} = hstore; 
    
    Vstore  = [Vstore Vbar]; 
    Vstore_no_mu = [Vstore_no_mu Vbar_no_mu];
    fprintf('Contact violation %.3f \n',norm(h));
%     addpoints(an,k_AL,norm(h));
%     drawnow
    hbar(:,k_AL) = h;
    if k_AL == 6 
        1==1
%         maxItr = 50; 
    end
    if params.DDP_OK
        norm_hbar(k_AL) = norm(h);
        if norm(h) < 1e-5 || k_AL > 4
            break; 
        end
        k_AL = k_AL + 1;
        mu1 = params.mu1 + params.mu2*h; % update lagrangian
        params.mu1 = mu1;
        params.mu2 = min(Beta_AL*params.mu2,1e8); % update penalty

        params.delta = max(params.delta*0.5,1e-7); 
        eps_list = 0.01* ones(1,20); 
%         eps_list(1) = 1; 
        params.eps = eps_list(k_AL);
        coeff_list = 0.7*ones(1,20)*10; 
        coeff_list(1)=3; coeff_list(2)=2; coeff_list(3)=1;
        params.coeff = coeff_list(k_AL - 1);
        umax_list = 34*ones(1,20); 
        umax_list(1)= 70; umax_list(2)=40; umax_list(3)=40;
        params.umax=umax_list(k_AL - 1);
        params.epsU = 600;

        
    else
        fprintf('DDP does not converge. Program terminates. \n');
        break;
    end
    %
%     if k_AL > 3
%         fprintf('=======================================\n');
%         fprintf('One Iteration Complete');
%         break
%     end
    
end
Time = toc(start);
Storage.Time = Time; 
Storage.rbtparams = robot_params; 
Storage.params = params; 
Storage.h0 = h0;

%%
%{
close all;
for idx = 1:k_AL
    fprintf('i: %i\n',idx);
    fprintf('Max Ubar: %.3f\n',max(max(Storage.ubar{idx})))
    fprintf('Min Ubar: %.3f\n',min(min( Storage.ubar{idx})))
end 

fname = strcat('ConstraintsOther_','iLQR','.gif');
params.filename = fname;
simulation(Storage.xbar{5}(1:params.q_size,:),robot_params,params,false);
fprintf('delta: %.3f\n',params.delta)
fprintf('mu1: %.3f\n',params.mu1)
fprintf('mu2: %.3f\n',params.mu2)

%
gg = 1:length(ubar); gg = gg*params.dt;
figure; hold on;
for idx=1:4
   plot(gg,Storage.ubar{5}(idx,:)); 
end 
grid on; grid minor
%
figure; hold on; plot(ybar(1,:,1)); plot(3*ybar(2,:,1),'DisplayName','0.7fz');legend;
figure; hold on; plot(ybar(1,:,2)); plot(3*ybar(2,:,2),'DisplayName','0.7fz');legend;
%

for i = 1:k_AL
    norm_hbar(i) = norm(hbar(:,i));
end
figure;
norm_hbar = [h0 norm_hbar(1:k_AL)];
plot(1:k_AL+1,norm_hbar,'linewidth',2);
xlabel('Outer loop iteration');
ylabel('$h(\mathbf{x})$','interpreter','latex');
set(gca,'fontsize',14);

%}

%% Create animation
%{
strMethod = switchMethod(Method);
strReg = switchRegularization(regularizationMethod);  

fname = strcat('figures/HeightAdded_','iLQR',num2str(iLQR),'method',strMethod,'Reg',strReg,'.gif');
params.filename = fname;
if strcmp(gait,"Trot") 
    simulation(xbar(1:params.q_size,:),robot_params,params);
else
    simulation_Jump(xbar(1:params.q_size,:),robot_params,params,true);
end
%}


function str = switchRegularization(regMethod)    
    switch regMethod 
        case 1 
            str = 'Traditional';
        case 2
            str = 'PointByPoint'; 
        case 3 
            str = 'TensPullBack';
        otherwise
            str = 'None';
    end
end

function str = switchMethod(method)    
    switch method 
        case 1 
            str = 'Explicit';
        case 2
            str = 'ExtMod'; 
        case 3 
            str = 'Tensor';
        otherwise
            str = 'None';
    end
end
end




