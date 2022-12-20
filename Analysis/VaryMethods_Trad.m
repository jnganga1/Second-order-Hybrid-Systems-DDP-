%% Run (full) iLQR and DDP (using Trad Regularization)
Itgr = 1; %DON'T CHANGE ME 
OtherTests.Run = false; 
OtherTests.randomInit = false;

%%
%Using iLQR method - no regularization
iLQR = 1; Method = 'none'; Reg =1;
iLQR_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);


%%
%using DDP, regularization: Trad, Method: ExtMod 
iLQR = 0; Method = 2; Reg = 1;
DDP_ExtMod_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);

%%
%using DDP, regularization: Trad, Method: Tensor 
iLQR = 0; Method = 3; Reg = 1;
DDP_Tens_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);
%%
save('MatData/VaryMethods_Trad_Data.mat');

%%
% close all; clearvars; clc; 
%%
%Making defaults makes life easier 
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)
%
%%
load('MatData/VaryMethods_Trad_Data.mat');


%Keep this order for clean code
colors = {'m','b','r'};
names = {'iLQR', 'mRNEAc DDP','Tensor DDP'};
Stuff = {iLQR_store,DDP_ExtMod_store,DDP_Tens_store}; 

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

%%

for i=1:length(Stuff)
    fprintf(['End Cost',names{i},':\n'])
%     vpa(Stuff{i}.Vbar{5}(end),1000)
    vpa(Stuff{i}.params.mu1,1000)
    vpa(Stuff{i}.params.mu2,1000)
%     Stuff{i}.params.delta
    
end
%%
figure;hold on; h ={};
for i=1:length(Stuff)
    a = [Stuff{i}.Vbar{5}];     
%     a = flipud(a);
    h{i}= semilogy(a- a(end),'DisplayName',names{i},'Color',colors{i}); 
end
legend

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
dd = 100*[1:8];
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
fprintf('Time: Tensor vs conventional FO %.2f\n',Time(3)/Time(1))
fprintf('Time: Tensor vs proposed %.2f\n',Time(3)/Time(2))
fprintf('Time: proposed vs conventional FO %.2f\n',Time(2)/Time(1))

fprintf('Prcnt Time: Tensor vs iLQR %.2f prcnt\n',100*(Time(3) - Time(1) )/ Time(1))
fprintf('Prcnt Time: Proposed vs iLQR %.2f prcnt\n',100*(Time(2) - Time(1) )/ Time(1))

bb = yy(:,10)
fprintf('Time: Tensor vs conventional FO %.2f\n',bb(3)/bb(1))
fprintf('Time: Tensor vs proposed %.2f\n',bb(3)/bb(2))
fprintf('Time: proposed vs conventional FO %.2f\n',bb(2)/bb(1))

disp('Time'); 
Time
%% Do some simulations and save the gif
for i =1:length(Stuff)
    fname = ['Data_BoundingNewNew/Results/',names{i},'iLQR','_Trad','.gif'];
    Stuff{i}.params.filename = fname;
    simulation_Jump(Stuff{i}.xbar{5},...
        Stuff{i}.ybar{5},Stuff{i}.rbtparams,Stuff{i}.params,false)
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
ybar =  Stuff{1}.ybar{5};
len = length(ybar);
tspan = dt*(0:len-1);



figure;
subplot(2,1,1)
plot(tspan,ybar(1,:,1),'DisplayName','$F_{x}$'); hold on; 
xlabel('time s','Interpreter','latex');
ylabel('GRF forces, $N$','interpreter','latex');
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
ylabel('GRF forces, $N$','interpreter','latex');
title('DDP: $F_{\rm{bc}}$','interpreter','latex'); legend;

%% Plot some ubar 
figure;
h={};
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
ylabel(['u',num2str(idx),' (Nm)'],'Interpreter','latex'); 
end
bb=legend([h{:}],'Location','best'); 
bb.Box='off'
%
iAL = 5;

fprintf('Max U: %.2f\n', max(max(Stuff{1}.ubar{iAL})));
fprintf('Max U: %.2f\n', min(min(Stuff{1}.ubar{iAL})));
%%
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
%% Plot xbar 

figure; G = {};
for idx = 1:length(Stuff)
    for i = 1:7
        subplot(2,4,i); hold on;
        G{1}=plot(Stuff{idx}.xbar{5}(i,:),'LineWidth',2);
        hold on; 
        grid on; grid minor;
        ylabel(['x',num2str(i)],'Interpreter','latex');
        xlabel('Time','Interpreter','latex');
    end
end
sgtitle('${q_i}$','Interpreter','latex'); 

figure; G = {};
for idx = 1:length(Stuff)
    for i = 1:7
        subplot(2,4,i); hold on;
        G{1}=plot(Stuff{idx}.xbar{5}(i+7,:),'LineWidth',2);
        hold on; 
        grid on; grid minor;
        ylabel(['x',num2str(i+7)],'Interpreter','latex');
        xlabel('Time','Interpreter','latex');
    end
end
sgtitle('${\dot{q}_i}$','Interpreter','latex'); 

%% Plot ubar


figure; G = {};
for idx = 1:length(Stuff)
    for i = 1:4
        subplot(2,2,i); hold on;
        G{1}=plot(Stuff{idx}.ubar{5}(i,:),'LineWidth',2);
        hold on; 
        grid on; grid minor;
        ylabel(['u',num2str(i)],'Interpreter','latex');
        xlabel('Time','Interpreter','latex');
    end
end
sgtitle('${u}$','Interpreter','latex'); 


%%
dt = 0.001;
ybar =  Stuff{1}.ybar{5};
len = length(ybar);
tspan = dt*(0:len-1);



figure;
subplot(2,1,1)
plot(tspan,ybar(1,:,1),'DisplayName','$F_{x}$'); hold on; 
xlabel('Time (s)','Interpreter','latex');
ylabel('GRF forces, $N$','interpreter','latex');
plot(tspan,0.7*ybar(2,:,1),'r','DisplayName','$0.7 F_{z}$');
plot(tspan,-0.7*ybar(2,:,1),'r--','DisplayName','$-0.7 F_{z}$');
% title('DDP: $F_{\rm{fr}}$','interpreter','latex'); hold off
bb=legend; 
bb.Location='best'; 
bb.Box='off';
grid on; grid minor;

subplot(2,1,2)
h={};
names = {'iLQR','DDP'};
for idx = 1:1
%     subplot(2,2,idx)
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
ylabel(['u',num2str(idx),' (Nm)'],'Interpreter','latex'); 
end
bb=legend([h{:}],'Location','best'); 
bb.Box='off'


xlim([time(1) time(end)]);

