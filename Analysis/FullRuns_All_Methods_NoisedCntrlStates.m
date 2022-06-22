%% Add Weiner noise to ubar - Run full iLQR/DDP - using Trad Regularization

%First use small noise and run DDP with fixed AL params for a few
%iterations --> will serve as basis for bigger run and short runs

close all; %clear; clc; 
Itgr = 1;
OtherTests.Run = true;  
OtherTests.randomInit = false;
OtherTests.maxItr =4000;
OtherTests.dt = 1e-3;

OtherTests.BeforeNoise = true; %fixing AL parameters and just running with it

OtherTests.varyU  = true; %care about varying u in this round

iLQR_Nsd = {}; ExtMod_Nsd = {};
Exp_Nsd = {}; Tens_Nsd = {};

%Control Noise
scalarNse = 10; %small noise
N = 317; dt = 1e-3; 
dW = sqrt(dt)*randn(4,N);   
W = cumsum(dW)*scalarNse; 
OtherTests.noise = W ;
OtherTests.scalarNse = scalarNse;

%State Noise
StateScalar = 0.1; %small noise


N = 317; dt = 1e-3; 
dW = sqrt(dt)*randn(14,N);   
W = cumsum(dW)*StateScalar; 
% W(8:end,:) = 0 *W(8:end,:);
OtherTests.StateNoise = W;
OtherTests.StateScalar = StateScalar;



iLQR = 1; Method = 'none'; Reg=  1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method;
OtherTests.RegMethod = 1;

iLQR_Nsd_nowOpt = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR
1==1;
%
iLQR = 0; Method = 2; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 

ExtMod_Nsd_nowOpt = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod

save('MatData/diditwork_again.mat')
fprintf('Returning from full run')
return

%%
data = load('MatData/diditwork_again.mat')
figure
semilogy(data.ExtMod_Nsd_nowOpt.Vbar{:} - data.ExtMod_Nsd_nowOpt.Vbar{:}(end),'DisplayName','DDP')
hold on; 
semilogy(data.iLQR_Nsd_nowOpt.Vbar{:} - data.iLQR_Nsd_nowOpt.Vbar{:}(end),'DisplayName','iLQR')
legend


fprintf('Orig DDP End Cost:  %.6f\n', vpa(iLQR_Nsd_nowOpt.Vbar{:}(end),6))
fprintf('Orig iLQR End Cost:  %.6f\n', vpa(iLQR_Nsd_nowOpt.Vbar{:}(end),6))
fprintf('iLQR End Cost: %.6f\n', vpa(iLQR_Nsd_nowOpt.Vbar{:}(end),6))
fprintf('DDP End Cost: %.6f\n', vpa(ExtMod_Nsd_nowOpt.Vbar{:}(end),6))


%%
OtherTests.BeforeNoise = false;

Itgr = 1;
OtherTests.Run = true; 
OtherTests.randomInit = false;
OtherTests.maxItr =5000;
OtherTests.dt = 1e-3;
OtherTests.varyU  = true; %care about varying u in this round



%Control Noise
scalarNse = 40; %More Noise 
N = 317; dt = 1e-3; 
dW = sqrt(dt)*randn(4,N);   
W = cumsum(dW)*scalarNse; 
OtherTests.noise = W ;
OtherTests.scalarNse = scalarNse;

%State Noise
StateScalar = 1.0; %More Noise


N = 317; dt = 1e-3; 
dW = sqrt(dt)*randn(14,N);   
W = cumsum(dW)*StateScalar; 
% W(8:end,:) = 0 *W(8:end,:);
OtherTests.StateNoise = W;
OtherTests.StateScalar = StateScalar;


iLQR = 1; Method = 'none'; Reg=  1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method;
OtherTests.RegMethod = Reg;
iLQR_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR


iLQR = 0; Method = 2; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 

ExtMod_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod

figure; 
semilogy(iLQR_Nsd.Vbar{:}-iLQR_Nsd.Vbar{:}(end),'DisplayName','iLQR'); 
hold on
semilogy(ExtMod_Nsd.Vbar{:}-ExtMod_Nsd.Vbar{:}(end),'DisplayName','DDP'); 
legend


%}
save('MatData/VaryU_Noise_FullRun_Fixed_ReModded4.mat');
return 
%%
close all; clear; clc; 
% load('crcMatData/VaryU_Noise_FullRun_DDPstart_Again_50Scalar.mat');
load('crcMatData_New/Latest/VaryU_Noise_FullRun_WithReg.mat');
%%

load('NonConvexData/Trad_Methods_fixed_nonconvex_afterRun.mat');

load('Data_BoundingNewNew/VaryU_Noise_FullRun_Fixed.mat');

load('crcMatData_New/JohnAgain_NoConst.mat');
%%

set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

%Keep this order for clean code
colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};
Stuff = {iLQR_Nsd_nowOpt,ExtMod_Nsd_nowOpt};%Exp_Nsd, Tens_Nsd};

% gg = data.ExtMod_Nsd_ot;
gg = [Stuff{2}.OrigOptim.ExtMod_Nsd_nowOpt]; %.Vbar{:}];

stateNoise =Stuff{1}.OtherTests.StateNoise;
controlNoise = Stuff{1}.OtherTests.noise;

%Plot previous optim states vs now optim 
%Here it is q 
figure; G = {};
for i = 1:7
    subplot(2,4,i); hold on;
    G{1}=plot(gg.xbar{5}(i,:),'LineWidth',2,'DisplayName','Optimized Before Noise'); 
    hold on; 
    G{2}=plot(Stuff{1}.xbar{1}(i,:),'LineWidth',2,'DisplayName','Optimized After Noise'); 
    grid on; grid minor;
    ylabel(['x',num2str(i)],'Interpreter','latex');
    xlabel('Time','Interpreter','latex');
end
legend
sgtitle('${q_i}$','Interpreter','latex'); 


%Plot previous optim states vs now optim 
%Here it is qdot
figure; G = {};
for i = 1:7
    subplot(2,4,i); hold on;
    G{1}=plot(gg.xbar{5}(i+7,:),'LineWidth',2,'DisplayName','Optimized Before Noise'); 
    hold on; 
    G{2}=plot(Stuff{1}.xbar{1}(i+7,:),'LineWidth',2,'DisplayName','Optimized After Noise'); 
    grid on; grid minor;
    ylabel(['x',num2str(i+7)],'Interpreter','latex');
    xlabel('Time','Interpreter','latex');
end
legend
sgtitle('$\dot{q_i}$','Interpreter','latex'); 
%% Plotp previous optim control vs now optim 

figure; G = {};
for i = 1:4
    subplot(2,2,i); hold on;
    G{1}=plot(gg.ubar{5}(i,:),'LineWidth',2,'DisplayName','Optimized Before Noise'); 
    hold on; 
    G{2}=plot(Stuff{1}.ubar{1}(i,:),'LineWidth',2,'DisplayName','Optimized After Noise'); 
    grid on; grid minor;
    ylabel(['u',num2str(i)],'Interpreter','latex');
    xlabel('Time','Interpreter','latex');
end
legend
% sgtitle('$\dot{q_i}$','Interpreter','latex');

%%
%Plot noised state to be optimized 
figure; 
subplot(2,2,1); G = {};
for i = 1:7
    subplot(2,4,i); hold on;
    bb = gg.xbar{5}(i,:);
    G{1}=plot(bb,'LineWidth',2,'DisplayName','x'); 
    hold on; 
    G{2}=plot(bb+stateNoise(i,1:length(bb)),'LineWidth',2,'DisplayName','x+Noise'); 
    grid on; grid minor;
    ylabel(['x',num2str(i)],'Interpreter','latex');
    xlabel('Time','Interpreter','latex');
    h=gca; h=h.Children;
    set(gca,'Children',[flipud(h)])
end
legend

% sgtitle('$\dot{q_i}$','Interpreter','latex');

%%

% figure; G = {};
subplot(2,2,2); G = {};
for i = 1:4
    subplot(2,2,i); hold on;
    bb=gg.ubar{5}(i,:);
    G{1}=plot(bb,'LineWidth',2,'DisplayName','u'); 
    hold on; 
    G{2}=plot(bb+controlNoise(i,1:length(bb)),'LineWidth',2,'DisplayName','u+Noise'); 
    grid on; grid minor;
    ylabel(['u',num2str(i)],'Interpreter','latex');
    xlabel('Time','Interpreter','latex');
%       h=gca; h=h.Children;
%     set(gca,'Children',[flipud(h)])
end
legend
% sgtitle('$\dot{q_i}$','Interpreter','latex');




%%
% Vbar_prev = gg.Vbar{:};
Vbar_iLQR = iLQR_Nsd.OrigOptim.iLQR_store.Vbar{5}(end);
Vbar_DDP = iLQR_Nsd.OrigOptim.DDP_ExtMod_store.Vbar{5}(end);
Vbar_end = min(Vbar_iLQR,Vbar_DDP);
%%
figure; h={};
gg=data.ExtMod_Nsd_ot.Vbar{1}(end);
%
for i=1:length(Stuff)
    
    hh=Stuff{i}.Vbar{:};
    
    h{i} = semilogy(hh - gg(end),'DisplayName',names{i});
    hold on;
end
legend([h{:}]);
grid on; grid minor; 
xlabel('Iterations','Interpreter','latex');
ylabel('Sub-optimality','Interpreter','latex'); 


%% Plot for Paper 

figure; 
gg = [Stuff{2}.OrigOptim.DDP_ExtMod_store]; %.Vbar{:}];


subplot(2,2,2); G = {};
for i = 4
%     subplot(2,2,i); hold on;
    bb=gg.ubar{5}(i,:);
    G{1}=plot(bb,'LineWidth',2,'DisplayName','u'); 
    hold on; 
    G{2}=plot(bb+controlNoise(i,1:length(bb)),'LineWidth',2,'DisplayName','u+Noise'); 
    grid on; grid minor;
    ylabel(['u',num2str(i)],'Interpreter','latex');
    xlabel('Time','Interpreter','latex');
      h=gca; h=h.Children;
    set(gca,'Children',[flipud(h)])
end
legend
% sgtitle('$\dot{q_i}$','Interpreter','latex');

% figure; 
subplot(2,2,4); G = {};
for i = 7
%     subplot(2,4,i); hold on;
    bb = gg.xbar{5}(i,:);
    G{1}=plot(bb,'LineWidth',2,'DisplayName','x'); 
    hold on; 
    G{2}=plot(bb+stateNoise(i,1:length(bb)),'LineWidth',2,'DisplayName','x+Noise'); 
    grid on; grid minor;
    ylabel(['x',num2str(i)],'Interpreter','latex');
    xlabel('Time','Interpreter','latex');
    h=gca; h=h.Children;
    set(gca,'Children',[flipud(h)])
end
legend


gg=data.ExtMod_Nsd_ot.Vbar{1}(end);

subplot(2,2,[1 3]); h={};
for i=1:length(Stuff)
    
    hh=Stuff{i}.Vbar{:};
    
    h{i} = semilogy(hh - gg(end),'DisplayName',names{i});
    hold on;
end
legend([h{:}]);
grid on; grid minor; 
xlabel('Iterations','Interpreter','latex');
ylabel('Sub-optimality','Interpreter','latex'); 


%%

for i=1:length(Stuff)
    hh = Stuff{i}.VnoMu{:}; hold on; 
    
    % Suboptimality
    %     semilogy(hh - hh(end),'DisplayName',names{i},'Color',colors{i}); 
    
    %Convergence Rate
    H = hh - gg(end); 
%     suLinear = abs(  H(2:end)./  H(1:end-1));
    Quad = abs( H(2:end) ./ H(1:end-1).^2);
% %     
%     suLinear = abs(  H(1:end-1)./  H(2:end));
%     Quad = abs( H(1:end-1) ./ H(2:end).^2);
%     
    
    figure(2); hold on;
%     semilogy(suLinear,'--','DisplayName',['Linear ',names{i}],'Color',colors{i},...
%         'Markersize',2);
%     figure(3); hold on;
    semilogy(Quad,'DisplayName',['Quad ',names{i}],'Color',colors{i});
    
    % 
%     fprintf(names{i});
%     max(max(Stuff{i}.ubar{:}))
%     min(min(Stuff{i}.ubar{:}))
    
end 
%
Conv = {'Super-linear Convergence Rate',...
    'Quadratic Convergence Rate'};
for i =2 
    figure(i)
    xlabel('Iterations','Interpreter','latex');
    ylabel(Conv{i-1},'Interpreter','latex');
    grid on; grid minor;
    bb = gca;
    bb.YScale='log';
    legend;
end
figure(2); 
title('Definition 2','Interpreter','latex');

%%
%Keep this order for clean code
colors = {'r','b','k','c'};
% Stuff = {iLQR_Nsd5,ExtMod_Nsd5};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};
for i =1:length(Stuff)
    fname = ['Results/crcMatData_New/',names{i},'iLQR','_LargeNoiseControl','.gif'];
    Stuff{i}.params.filename = fname;
    simulation_Jump(Stuff{i}.xbar{1},...
        Stuff{i}.ybar{1},Stuff{i}.rbtparams,Stuff{i}.params,false)
%     norm(Stuff{i}.ybar{idx}(:))
end

