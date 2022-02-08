%% Run (full) iLQR and DDP (using TradMod Regularization)
Itgr = 1; %DON'T CHANGE ME 
OtherTests.Run = false; 
OtherTests.randomInit = false;
%
%Using iLQR method - no regularization
iLQR = 1; Method = 'none';
iLQR_store  = DDP_2DQuadruped(iLQR,Method,0,Itgr,OtherTests);
%
Reg = 4; 
%using DDP, regularization: TradMod, Method: ExtMod 
iLQR = 0; Method = 2;
DDP_ExtMod_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);
%
%using DDP, regularization: TradMod, Method: Explicit 
iLQR = 0; Method = 1; 
DDP_Exp_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);

%using DDP, regularization: TradMod, Method: ExtMod 
iLQR = 0; Method = 3;
DDP_Tens_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);


save('MatData/TradMod_Methods_compare.mat');

%%
close all; clear; clc; 
load('MatData/TradMod_Methods_compare.mat');

%Making defaults makes life easier 
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

load('TradMod_Methods_compare.mat');
%Keep this order for clean code
colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};
Stuff = {iLQR_store, DDP_ExtMod_store,DDP_Exp_store,DDP_Tens_store}; 

% Cost vs Iterations
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
maxItr = 100;
for i=1:length(Stuff)
    for idx =1:length(Stuff{i}.Vbar)
            ss = (idx-1)*maxItr; 
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
    fname = ['Results/',names{i},'iLQR','_TradMod','.gif'];
    Stuff{i}.params.filename = fname;
    simulation_Jump(Stuff{i}.xbar{5},...
        Stuff{i}.ybar{5},Stuff{i}.rbtparams,Stuff{i}.params,true)
%     norm(Stuff{i}.ybar{idx}(:))
end

%% Plot some ybar - in x and z directions
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
%% Plot Ybar - Friction cone wise

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

hold all
x = [1 1]*maxItr; 
h = gca; hold on 
for i= 1:length(ilqr_Vstore)/100
    l=plot(i*x,h.YLim,'k','LineWidth',2); 
end
hold off
h = legend;
h.Location = 'best';
h.String(5:end) = '';
