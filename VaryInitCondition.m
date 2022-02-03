%% Run Tests Based on Limited iterations of DDP  - makes figure B,C,D
%Using iLQR method - no regularization
iLQR = 1; Method = 'none';

OtherTests.randomInit = false;

OtherTests.Run = true; 
dt_lst = logspace(-5,-0.01,10);
OtherTests.maxItr = 10;

ExpEurBig = {}; ExpMidBig ={};
ImpEulrBig = {}; ImpMidBig = {}; 
for i = 1:length(dt_lst)
    OtherTests.dt = dt_lst(i); 
    ExpEurBig{i} = DDP_2DQuadruped(iLQR,Method,0,1,OtherTests);
    ImpEulrBig{i} = DDP_2DQuadruped(iLQR,Method,0,2,OtherTests);
    ExpMidBig{i} = DDP_2DQuadruped(iLQR,Method,0,3,OtherTests);
    ImpMidBig{i} = DDP_2DQuadruped(iLQR,Method,0,4,OtherTests);
    close all;
end 


%% Random Initialization 
%for loop here with randomization
OtherTests.randomInit = true;
OtherTests.Run = true; 
OtherTests.maxItr = 10;
OtherTests.dt = 1e-3;
iRepts = 50; 


rndValue = [];
ExpEur_Rndm = {}; ExpMid_Rndm ={};
ImpEulr_Rndm = {};ImpMid_Rndm = {}; 
a=-1.5;b=1.5; 
for i = 1:iRepts     
    r = a + (b-a).*rand(1);
    rndValue(end+1) = r;
    OtherTests.rndValue = r;
    ExpEur_Rndm{i} = DDP_2DQuadruped(iLQR,Method,0,1,OtherTests);
    ImpEulr_Rndm{i} = DDP_2DQuadruped(iLQR,Method,0,2,OtherTests);
    ExpMid_Rndm{i} = DDP_2DQuadruped(iLQR,Method,0,3,OtherTests);
    ImpMid_Rndm{i} = DDP_2DQuadruped(iLQR,Method,0,4,OtherTests);
    close all;
end 

% save('VaryInitCond');
%% 
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)
%

%Keep this order for clean code
colors = {'r','b','k','c'};
names = {'Explicit Euler', 'Implicit Euler',...
    'Explicit Midpoint', 'Implicit Midpoint'};  
markers = {'s','*','+','o'};

%% Plot EndCost Vs dt

Exp_EndCost = []; Imp_EndCost = [];
ImpMid_EndCost  = []; ExpMid_EndCost= [];

Exp_Time = []; Imp_Time = [];
ImpMid_Time  = []; ExpMid_Time= [];

for i = 1:length(dt_lst)
    Exp_EndCost(end+1) = ExpEurBig{i}.Vbar{:}(end);
    Imp_EndCost(end+1) = ImpEulrBig{i}.Vbar{:}(end);
    ExpMid_EndCost(end+1) = ExpMidBig{i}.Vbar{:}(end);
    ImpMid_EndCost(end+1) = ImpMidBig{i}.Vbar{:}(end);
    
    Exp_Time(end+1) = ExpEurBig{i}.Time;
    Imp_Time(end+1) = ImpEulrBig{i}.Time; 
    ExpMid_Time(end+1) = ExpMidBig{i}.Time;
    ImpMid_Time(end+1) = ImpMidBig{i}.Time;
end

Stuffed = [Exp_EndCost;Imp_EndCost;ExpMid_EndCost;ImpMid_EndCost];

figure; hold on;
for idx = 1:size(Stuffed,1) 
    plot(dt_lst,Stuffed(idx,:),'DisplayName',names{idx},...
        'color',colors{idx});%,'Marker',markers{idx});
end
grid on; grid minor;
hh = legend;
hh.Location = 'best';
xlabel('dt (s)','Interpreter','latex'); 
ylabel('Cost at end of 10 Iterations','Interpreter','latex');
bb=gca; 
bb.XLim(2)=.0017;

%% Plot Time vs dt 

Stuffed_Time = [Exp_Time;Imp_Time;ExpMid_Time;ImpMid_Time];
ind = ~isnan(Stuffed);

figure; hold on;
for idx = 1:size(Stuffed_Time,1) 
    a = Stuffed_Time(idx,:); 
    a = a(ind(idx,:));
    plot(dt_lst(ind(idx,:)),a,'DisplayName',names{idx},...
        'color',colors{idx});%,'Marker',markers{idx});
end
grid on; grid minor; 
hh = legend;
hh.Location = 'best';
xlabel('dt (s)','Interpreter','latex'); 
ylabel('Time(s) of 10 Iterations','Interpreter','latex');
%
bb=gca; 
bb.XLim(2)=.0017;

%% Plot Success vs dt 

ind = isnan(Stuffed);
Succ = Stuffed; 
Succ(~ind) = 1; 
Succ(ind)  = 0;

figure; hold on
for idx = 1:size(Succ,1) 
    plot(dt_lst,Succ(idx,:),'DisplayName',names{idx},...
        'color',colors{idx},'LineStyle','none','Marker',markers{idx});
end
%
bb=gca;
bb.XLim(2) = 0.02;grid on; grid minor; 
yticks([0 1])
yticklabels({'Fail','Success'});
xlabel('dt (s)','Interpreter','latex'); 
legend

%% Randomization of init condition 

Exp_Rnd= []; Imp_Rnd = []; ExpMid_Rnd=[]; ImpMid_Rnd=[];
for i = 1:50
%     Exp_Rnd(end+1) = ExpEur_Rndm{i}.Vbar{:}(end); 
%     Imp_Rnd(end+1) = ImpEulr_Rndm{i}.Vbar{:}(end);
%     ExpMid_Rnd(end+1) = ExpMid_Rndm{i}.Vbar{:}(end);
%     ImpMid_Rnd(end+1) = ImpMid_Rndm{i}.Vbar{:}(end);
    
    Exp_Rnd(end+1) = norm(ExpEur_Rndm{i}.ybar{:}(:));
    Imp_Rnd(end+1) = norm(ImpEulr_Rndm{i}.ybar{:}(:));
    ExpMid_Rnd(end+1) = norm(ExpMid_Rndm{i}.ybar{:}(:));
    ImpMid_Rnd(end+1) = norm(ImpMid_Rndm{i}.ybar{:}(:));
    
end

Stuffed_Rnd = [Exp_Rnd; Imp_Rnd; ExpMid_Rnd; ImpMid_Rnd];
figure; hold on; 
for idx=1:size(Stuffed_Rnd,1)
    plot(rndValue,Stuffed_Rnd(idx,:),'DisplayName',names{idx},...
        'Marker',markers{idx},'MarkerSize',3,'LineStyle','none',...
        'color',colors{idx}); 
end
grid on; grid minor; 
hh = legend;
hh.Location = 'best';
xlabel('Variation to $\dot{v}_{Trunk} \; += 0.5$','Interpreter','latex')
ylabel('Cost at end of 10 Iterations','Interpreter','latex');
