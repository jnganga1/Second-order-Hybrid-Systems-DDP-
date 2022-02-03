%% Run (full) iLQR and DDP (using Trad Regularization)
Itgr = 1; %DON'T CHANGE ME 
OtherTests.Run = false; 
OtherTests.randomInit = false;

%Using iLQR method - no regularization
iLQR = 1; Method = 'none';
iLQR_store  = DDP_2DQuadruped(iLQR,Method,0,Itgr,OtherTests);

%using DDP, regularization: Trad, Method: ExtMod 
iLQR = 0; Method = 2; Reg = 1;
DDP_ExtMod_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);

%using DDP, regularization: Trad, Method: Explicit 
iLQR = 0; Method = 1; Reg = 1;
DDP_Exp_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);

%using DDP, regularization: Trad, Method: ExtMod 
iLQR = 0; Method = 3; Reg = 1;
DDP_Tens_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);


save('Trad_Methods_compare.mat');

%% Add noise to ubar - Run iLQR/DDP for 20 Iterations 
close all; clear; clc; 
% Add noise to ubar -- some weiner noise.
% Run iLQR/DDP for 50 Iterations and get time - makes figure A 
Itgr = 1;
OtherTests.Run = true; 
OtherTests.randomInit = false;
OtherTests.maxItr = 1000;
OtherTests.dt = 1e-3;

OtherTests.varyU  = true; %care about varying u in this round

iLQR_Box = {}; ExtMod_Box ={};
Exp_Box = {}; Tens_Box= {};

% scaler = [-3 -1 1 2 3];
% scaler = [-2:2];
scaler = linspace(50,50,7);
% scaler = [-2 0 1 2];
iRepts = length(scaler);
% 
% figure;
% for i = 1:iRepts    
% %     N = length(ubar);
%     N = 585; dt = 1e-3; 
%     dW = sqrt(dt)*randn(4,N);   
%     W = cumsum(dW) *  scaler(i);
%     fprintf('scaler: %i, maxNoise: %f\n',scaler(i),max(W(1,:)));
%     plot(W(1,:),'DisplayName',num2str(scaler(i))); hold on
% end
% legend

cWork = {};

iScal = 7;
for idx = 1:iScal
    for i = 1:iRepts    
    %     N = length(ubar);
%%
        N = 585; dt = 1e-3; 
        dW = sqrt(dt)*randn(4,N);   
        W = cumsum(dW) *  scaler(i);

        OtherTests.noise = W; 
        cWork{idx,i} = {W,scaler(i)}; 

   
        iLQR = 1; Method = 'none'; Reg=  'none';
        OtherTests.iLQR = iLQR; OtherTests.Method = Method;
        OtherTests.RegMethod = 'none';

        iLQR_Box{idx,i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR

        iLQR = 0; Method = 2; Reg = 1;
        OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
        OtherTests.RegMethod = Reg; 

        ExtMod_Box{idx,i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod
        
  %%      
        figure; 
        hh=iLQR_Box{1}.VnoMu{:};
        semilogy(hh -hh(end),'DisplayName','iLQR'); hold on;
        gg=ExtMod_Box{1}.VnoMu{:};
        semilogy(gg - gg(end),'DisplayName','DDP');
        legend
        figure; 
        plot(hh,'DisplayName','iLQR'); hold on;
        plot(gg,'DisplayName','DDP');
        
        %
        figure; 
        gg2 = gg - gg(end);
        hh2 = hh - hh(end);
        subplot(2,1,1);
        H = abs( sqrt(hh2(1:end-1)) ./ hh2(2:end));
        semilogy(H,'-o','DisplayName','iLQR quadratic'); hold on; 
        G = abs( sqrt(gg2(1:end-1)) ./ gg2(2:end));
        semilogy(G,'-o','DisplayName','DDP quadratic');
        legend; grid on; grid minor;
        
        subplot(2,1,2);
        H = abs( hh2(1:end-1)) ./ hh2(2:end);
        semilogy(H,'-o','DisplayName','iLQR_linear'); hold on; 
        G = abs(gg2(1:end-1) ./ gg2(2:end));
        semilogy(G,'-o','DisplayName','DDP_linear');
        legend; grid on; grid minor;
        
        
        data= load('Trad_Methods_compare.mat');
        fprintf('Optimal basin:\n');
        Vend = data.iLQR_store.Vbar{5}(end)
        fprintf('Achieved iLQR optimal basin:\n');
        hh(end)        
        fprintf('Achieved iLQR optimal basin:\n');
        gg(end)
        
        
        
%%
        iLQR = 0; Method = 1; Reg = 1;
        OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
        OtherTests.RegMethod = Reg; 

        Exp_Box{idx,i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Exp
        iLQR = 0; Method = 3; Reg = 1;

        OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
        OtherTests.RegMethod = Reg; 
        Tens_Box{idx,i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Tens
    end 
end

figure; 
plot(iLQR_Box{1}.Vbar{:},'DisplayName','iLQR'); hold on; 
plot(ExtMod_Box{1}.Vbar{:},'DisplayName','ExtMod');
plot(Exp_Box{1}.Vbar{:},'DisplayName','Exp');
plot(Tens_Box{1}.Vbar{:},'DisplayName','Tens');
legend

save('VaryU_Noise_MultiRuns.mat');
%% Add Weiner noise to ubar - Run full iLQR/DDP

close all; clear; clc; 
Itgr = 1;
OtherTests.Run = true; 
OtherTests.randomInit = false;
OtherTests.maxItr = 150;
OtherTests.dt = 1e-3;

OtherTests.varyU  = true; %care about varying u in this round

iLQR_Nsd = {}; ExtMod_Nsd = {};
Exp_Nsd = {}; Tens_Nsd = {};

N = 585; dt = 1e-3; 
dW = sqrt(dt)*randn(4,N);   
W = cumsum(dW); 

OtherTests.noise = W; 


iLQR = 1; Method = 'none'; Reg=  'none';
OtherTests.iLQR = iLQR; OtherTests.Method = Method;
OtherTests.RegMethod = 'none';

iLQR_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR

iLQR = 0; Method = 2; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 

ExtMod_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod

iLQR = 0; Method = 1; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 
Exp_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Exp

iLQR = 0; Method = 3; Reg = 1;
OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
OtherTests.RegMethod = Reg; 
Tens_Nsd = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Tens

save('VaryU_Noise_FullRun_DDPstart.mat');
%%
figure;
hh=iLQR_Nsd.Vbar{:};
semilogy(hh - hh(end),'DisplayName','iLQR')
hold on;
gg=ExtMod_Nsd.Vbar{:};
semilogy(gg - gg(end),'DisplayName','DDP'); hold off
legend

%%  Vary Initial Condition and Run for 20 Iterations
close all; clear; clc; 
% Add noise to ubar -- some weiner noise.
% Run iLQR/DDP for 50 Iterations and get time - makes figure A 
Itgr = 1;
OtherTests.Run = true; 
% OtherTests.randomInit = true;
OtherTests.maxItr = 20;
OtherTests.dt = 1e-3;

OtherTests.varyU  = false; %Dont care about varying u in this round


iLQR_dfVd = {}; ExtMod_dfVd  ={};
Exp_dfVd = {}; Tens_dfVd = {};

% scaler = [-3 -1 1 2 3];
% scaler = [-2:2];
vd_list = linspace(0,3.2,7);
iRepts = length(vd_list);


for i = 1:iRepts    
    OtherTests.vd = vd_list(i);
    
    iLQR = 1; Method = 'none'; Reg=  'none';
    OtherTests.iLQR = iLQR; OtherTests.Method = Method;
    OtherTests.RegMethod = 'none';
    
    iLQR_dfVd{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR
    
    iLQR = 0; Method = 2; Reg = 1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 
    
    ExtMod_dfVd{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod
    
    iLQR = 0; Method = 1; Reg = 1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 
    
    Exp_dfVd{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Exp
    
    iLQR = 0; Method = 3; Reg = 1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 
    Tens_dfVd{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Tens
end 

save('VaryVd.mat');


%%

%simulation(ExpEur.xbar{2},ExpEur.rbtparams,ExpEur.params,false)

set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)
%
%

load('Trad_Methods_compare.mat');
%Keep this order for clean code
colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};
Stuff = {iLQR_store, DDP_ExtMod_store,DDP_Exp_store,DDP_Tens_store}; 

%% Cost vs Iterations
figure;hold on; h ={};
for i=1:length(Stuff)
    a = [Stuff{i}.Vbar{1:5}];     
    h{i}= plot(a,'DisplayName',names{i},'Color',colors{i}); 
end
grid on; grid minor;legend; 
G=gca; G.YScale = 'log';
xlabel('Iterations','Interpreter','latex'); 
ylabel('Cost','Interpreter','latex'); 

%% Cost at each update vs Iterations
markers = {'s','*','+','o'}; 
figure;hold on; h={}; 
for i=1:length(Stuff)
    for idx =1:length(Stuff{i}.Vbar)
            ss = (idx-1)*100; 
            a =[Stuff{i}.Vbar{idx}];
            a = a - a(end); 
            xx = ss:ss+length(a)-1;        
            h{i}= plot(xx,a,'DisplayName',names{i},'Color',colors{i});
            if idx > 3 
                h{i}= plot(xx,a,'DisplayName',names{i},'Color',colors{i},...
                    'LineWidth',1,'Marker',markers{i}); 
            end
    end
end
grid on; grid minor;
G=gca; G.YScale = 'log';
xlabel('Iterations','Interpreter','latex'); 
ylabel('Cost','Interpreter','latex');
dd = 100*[1:4];
for i=1:length(dd)
    ssLine = dd(i)*ones(1,10); 
    bb = gca; bb = linspace(bb.YLim(1),bb.YLim(2),10);
    plot(ssLine,bb,'m--')
end
legend([h{:}],'Location','best'); 
%% Number of iterations per update and Total Time
Iters = []; Time = zeros(length(Stuff),1);

CostEnd=[];
for i=1:length(Stuff)
    for idx= 1:length(Stuff{i}.Vbar)-2
        Iters(i,idx) = length(Stuff{i}.Vbar{idx});
        CostEnd(i,idx) = Stuff{i}.Vbar{idx}(end); 
    end
    Time(i)=Stuff{i}.Time;
end
for i=1:length(Stuff) 
    fprintf('Iterations of %s',names{i}); 
    CostEnd(i,:)
end
disp('Time'); 
Time
%% Do some simulations and save the gif
for i =1:length(Stuff)
    fname = ['Results/',names{i},'iLQR','.gif'];
    Stuff{i}.params.filename = fname;
    simulation_Jump(Stuff{i}.xbar{5},...
        Stuff{i}.ybar{5},Stuff{i}.rbtparams,Stuff{i}.params,true)
%     norm(Stuff{i}.ybar{idx}(:))
end

%% Plot some ybar 
figure; hold on
subplot(2,1,1);
for i =1:length(Stuff) 
    for idx=1
        fz = Stuff{i}.ybar{idx}(idx,:,1);
        dt = 1e-3; time = 0:dt:(length(fz)-1)*dt;
        h{i} =plot(time,fz,...
            'DisplayName',names{i},'Color',colors{i});hold on;
        grid on; grid minor;
    end 
end
ylabel('Force (Nm)','Interpreter','latex'); 
% legend([h{:}],'Location','best'); 
subplot(2,1,2)
for i =1:length(Stuff) 
    for idx=2
        fz = Stuff{i}.ybar{idx}(idx,:,1);
        dt = 1e-3; time = 0:dt:(length(fz)-1)*dt;
        h{i} =plot(time,fz,...
            'DisplayName',names{i},'Color',colors{i});hold on;
        grid on; grid minor;
    end 
end
legend([h{:}],'Location','best'); 
xlabel('Time(s)','Interpreter','latex');
ylabel('Force (Nm)','Interpreter','latex'); 
%% Plot Ybar 

dt = 0.001;
len = 585;
tspan = dt*(0:len-1);

ybar =  Stuff{2}.ybar{5};
figure;
subplot(2,1,1)
plot(tspan,ybar(1,:,1),'DisplayName','$F_{x}$'); hold on; 
xlabel('time s','Interpreter','latex');
ylabel('$GRF force$','interpreter','latex');
plot(tspan,0.7*ybar(2,:,1),'r','DisplayName','$0.7 F_{z}$');
plot(tspan,-0.7*ybar(2,:,1),'r--','DisplayName','$-0.7 F_{z}$');
title('DDP: $F_{\rm{fr}}$','interpreter','latex'); hold off
% hl = legend('show');
% set(hl, 'Interpreter','latex')

subplot(2,1,2)
plot(tspan,ybar(1,:,2),'DisplayName','$F_{x}$'); hold on;
xlabel('time s','interpreter','latex');
title('$F_{\rm{bc}}$','interpreter','latex');
plot(tspan,0.7*ybar(2,:,2),'r','DisplayName','$0.7 F_{z}$');
plot(tspan,-0.7*ybar(2,:,2),'r--','DisplayName','$-0.7 F_{z}$');
xlabel('time s','interpreter','latex');
ylabel('$GRF forces$','interpreter','latex');
title('DDP: $F_{\rm{bc}}$','interpreter','latex'); legend;

    
   

%% Plot some ubar 
figure;
for idx = 1:4
    subplot(2,2,idx)
    for i =1:length(Stuff) 
        for id=idx
            u1 = Stuff{i}.ubar{5}(idx,:);
            dt = 1e-3; time = 0:dt:(length(u1)-1)*dt;
            h{i} =plot(time,u1,...
                'DisplayName',names{i},'Color',colors{i});hold on%,...
    %             'Marker',none);hold on;
        end 
    end
    bb=gca;
    maxu = 34;
    plot([bb.XLim(1) bb.XLim(2)],[-maxu -maxu],'k--');
    plot([bb.XLim(1) bb.XLim(2)],[maxu maxu],'k--');
    

grid on; 
grid minor;
xlabel('Time(s)','Interpreter','latex');
ylabel(['u',num2str(idx),' (N)'],'Interpreter','latex'); 
end
bb=legend([h{:}],'Location','best'); 
bb.Box='off'

% bb=gca;
% bb.XLim(1)=0.07; bb.XLim(2)=0.15;
%% Plot box plot of time
% close all; clear; clc;
load('VaryU_Noise_part2.mat');


colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};

iLQR_Box_Time = [];  ExtMod_Box_Time = [];
Exp_Box_Time  = [];  Tens_Box_Time = [];

ExtMod_normU = []; iLQR_normU = [];
Exp_normU = []; Tens_normU = [];
for i=1:iRepts
    iLQR_Box_Time(end+1) = iLQR_Box{i}.Time;     
    ExtMod_Box_Time(end+1) =  ExtMod_Box{i}.Time; 
    Exp_Box_Time(end+1)  = Exp_Box{i}.Time;     
    Tens_Box_Time(end+1) = Tens_Box{i}.Time;
    
    iLQR_normU(end+1) = norm(iLQR_Box{i}.ubar{:});
    ExtMod_normU(end+1) = norm(ExtMod_Box{i}.ubar{:});
    Exp_normU(end+1) = norm(Exp_Box{i}.ubar{:});
    Tens_normU(end+1) = norm(Tens_Box{i}.ubar{:});
   
end

Boxed_Time = [iLQR_Box_Time;ExtMod_Box_Time;...
    Exp_Box_Time;Tens_Box_Time]';
figure; 
boxplot(Boxed_Time,names); 
grid on; grid minor
ylabel('Time(s), 20 Iterations','Interpreter','latex');
bb=gca;
bb.XAxis.TickLabelInterpreter='latex'; 
bb.XTickLabelRotation=15;
bb.YScale = 'log';


%% Plot box plot of EndCost 


iLQR_Box_EndCost = [];  ExtMod_Box_EndCost = [];
Exp_Box_EndCost  = [];  Tens_Box_EndCost = [];

iLQR_length =[]; ExtMod_length =[];
Exp_length =[]; Tens_length = [];
for i=1:iRepts
    iLQR_Box_EndCost(end+1) = iLQR_Box{i}.Vbar{:}(end);     
    ExtMod_Box_EndCost(end+1) =  ExtMod_Box{i}.Vbar{:}(end); 
    Exp_Box_EndCost(end+1)  = Exp_Box{i}.Vbar{:}(end);     
    Tens_Box_EndCost(end+1) = Tens_Box{i}.Vbar{:}(end);
        
    iLQR_length(end+1) = length(iLQR_Box{i}.Vbar{:});   
    ExtMod_length(end+1) =  length(ExtMod_Box{i}.Vbar{:}); 
    Exp_length(end+1)  = length(Exp_Box{i}.Vbar{:});    
    Tens_length(end+1) = length(Tens_Box{i}.Vbar{:});
   
end


Boxed_EndCost = [iLQR_Box_EndCost;ExtMod_Box_EndCost;...
    Exp_Box_EndCost;Tens_Box_EndCost]';

iLQR_Box_EndCost(find(iLQR_length < 5)) = nan;
ExtMod_Box_EndCost(find(iLQR_length < 5)) = nan;
Exp_Box_EndCost(find(iLQR_length < 5)) = nan;
Tens_Box_EndCost(find(iLQR_length < 5)) = nan;

figure; 
boxplot(Boxed_EndCost,names); 
grid on; grid minor
ylabel('EndCost, 20 Iterations','Interpreter','latex');
bb=gca;
bb.XAxis.TickLabelInterpreter='latex'; 
bb.XTickLabelRotation=15;

%%
figure; hold on;
for i=1:4
    plot(scaler,Boxed_EndCost(:,i),'o','DisplayName',names{i});  
end
grid on; grid minor
ylabel('EndCost, 20 Iterations','Interpreter','latex');
xlabel('Noise Scaler');
bb=gca;
bb.XAxis.TickLabelInterpreter='latex'; 
legend
bb.YScale = 'log'

%%
figure; 
plot(scaler,iLQR_normU,'o','DisplayName','iLQR'); hold on
plot(scaler,ExtMod_normU,'o','DisplayName','ExtMod');
plot(scaler,Exp_normU,'o','DisplayName','Exp'); 
plot(scaler,Tens_normU,'o','DisplayName','Tens');
legend

%% 
close all; clear; clc; 
load('VaryU_Noise_FullRun_DDPstart.mat');

set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

%Keep this order for clean code
colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};
Stuff = {iLQR_Nsd, ExtMod_Nsd, Exp_Nsd, Tens_Nsd};

figure;
for i=1:length(Stuff)
    hh = Stuff{i}.VnoMu{:}; hold on; 
    
    % Suboptimality
    %     semilogy(hh - hh(end),'DisplayName',names{i},'Color',colors{i}); 
    
    %Convergence Rate
    H = hh - hh(end); 
    H = abs( H(1:end-1) ./ H(2:end));
    semilogy(H,'DisplayName',names{i},'Color',colors{i}); 
    
    % 
%     fprintf(names{i});
%     max(max(Stuff{i}.ubar{:}))
%     min(min(Stuff{i}.ubar{:}))
    
end 
xlabel('Iterations','Interpreter','latex');
ylabel('Convergence','Interpreter','latex');
grid on; grid minor;
bb = gca;
bb.YScale = 'log'; 
legend;

%%

load('VaryVd.mat'); 


colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};

iLQR_Box_Time = [];  ExtMod_Box_Time = [];
Exp_Box_Time  = [];  Tens_Box_Time = [];


for i=1:iRepts
    iLQR_Box_Time(end+1) = iLQR_dfVd{i}.Time;     
    ExtMod_Box_Time(end+1) =  ExtMod_dfVd{i}.Time; 
    Exp_Box_Time(end+1)  = Exp_dfVd{i}.Time;     
    Tens_Box_Time(end+1) = Tens_dfVd{i}.Time;
 
end

Boxed_Time = [iLQR_Box_Time;ExtMod_Box_Time;...
    Exp_Box_Time;Tens_Box_Time]';
figure; 
boxplot(Boxed_Time,names); 
grid on; grid minor
ylabel('Time(s), 20 Iterations','Interpreter','latex');
bb=gca;
bb.XAxis.TickLabelInterpreter='latex'; 
bb.XTickLabelRotation=15;
bb.YScale = 'log';


