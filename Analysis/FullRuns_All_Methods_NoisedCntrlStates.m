%% Add Weiner noise to ubar - Run full iLQR/DDP - using Trad Regularization
close all; clear; clc; 
Itgr = 1;
OtherTests.Run = true; 
OtherTests.randomInit = false;
OtherTests.maxItr = 1000;
OtherTests.dt = 1e-3;

OtherTests.varyU  = true; %care about varying u in this round

iLQR_Nsd = {}; ExtMod_Nsd = {};
Exp_Nsd = {}; Tens_Nsd = {};

%Control Noise
scalarNse = 200;
N = 585; dt = 1e-3; 
dW = sqrt(dt)*randn(4,N);   
W = cumsum(dW)*scalarNse; 
OtherTests.noise = W;
OtherTests.scalarNse = scalarNse;

%State Noise
StateScalar = 0;
N = 585; dt = 1e-3; 
dW = sqrt(dt)*randn(14,N);   
W = cumsum(dW)*StateScalar; 
W(8:end,:) = 0 *W(8:end,:);
OtherTests.StateNoise = W;
OtherTests.StateScalar = StateScalar;




iLQR = 1; Method = 'none'; Reg=  'none';
OtherTests.iLQR = iLQR; OtherTests.Method = Method;
OtherTests.RegMethod = 'none';

iLQR_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR
%%
iLQR = 0; Method = 2; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 

ExtMod_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod
%
iLQR = 0; Method = 1; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 
Exp_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Exp

iLQR = 0; Method = 3; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 
Tens_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Tens
%}
save('MatData/VaryU_Noise_FullRun_DDPstart_Again.mat');

%%
close all; clear; clc; 
load('MatData/VaryU_Noise_FullRun_DDPstart_Again.mat');

set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

%Keep this order for clean code
colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};
Stuff = {iLQR_Nsd, ExtMod_Nsd};% Exp_Nsd, Tens_Nsd};


figure(1); h={};
for i=1:length(Stuff) 
    hh=Stuff{i}.Vbar{:};
    h{i} = semilogy(hh - hh(end),'DisplayName',names{i});
    hold on;
end
legend([h{:}]);
grid on; grid minor; 
xlabel('Iterations','Interpreter','latex');
ylabel('Sub-optimality','Interpreter','latex'); 
%

for i=1:length(Stuff)
    hh = Stuff{i}.VnoMu{:}; hold on; 
    
    % Suboptimality
    %     semilogy(hh - hh(end),'DisplayName',names{i},'Color',colors{i}); 
    
    %Convergence Rate
    H = hh - hh(end); 
    suLinear = abs( H(1:end-1) ./ H(2:end));
    Quad = abs( H(1:end-1) ./ H(2:end).^2);
    figure(2); hold on;
    semilogy(suLinear,'DisplayName',names{i},'Color',colors{i});
    figure(3); hold on;
    semilogy(Quad,'DisplayName',names{i},'Color',colors{i});
    
    % 
%     fprintf(names{i});
%     max(max(Stuff{i}.ubar{:}))
%     min(min(Stuff{i}.ubar{:}))
    
end 
%
Conv = {'Super-linear Convergence Rate',...
    'Quadratic Convergence Rate'};
for i =2:3 
    figure(i)
    xlabel('Iterations','Interpreter','latex');
    ylabel(Conv{i-1},'Interpreter','latex');
    grid on; grid minor;
    bb = gca;
    bb.YScale='log';
    legend;
end

