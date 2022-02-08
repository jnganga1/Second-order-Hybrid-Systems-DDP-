%% Add noise to ubar - Run iLQR/DDP for 20 Iterations 
close all; clear; clc; 
% Add noise to ubar -- some weiner noise.
% Run iLQR/DDP for 50 Iterations and get time - makes figure A 
Itgr = 1;
OtherTests.Run = true; 
OtherTests.randomInit = false;
OtherTests.maxItr = 20;
OtherTests.dt = 1e-3;

OtherTests.varyU  = true; %care about varying u in this round

iLQR_Box = {}; ExtMod_Box ={};
Exp_Box = {}; Tens_Box= {};

% scaler = [-3 -1 1 2 3];
% scaler = [-2:2];
scaler = linspace(0,200,15);
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

iScal = 3;
for idx = 1:iScal
    for i = 1:iRepts    
    %     N = length(ubar);
%
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

% figure; 
% plot(iLQR_Box{1}.Vbar{:},'DisplayName','iLQR'); hold on; 
% plot(ExtMod_Box{1}.Vbar{:},'DisplayName','ExtMod');
% plot(Exp_Box{1}.Vbar{:},'DisplayName','Exp');
% plot(Tens_Box{1}.Vbar{:},'DisplayName','Tens');
% legend

save('MatData/VaryU_Noise_MultiRuns_Again.mat');

%%
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

load('MatData/VaryU_Noise_MultiRuns_Again.mat');


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
%% Plot box plot of time
% close all; clear; clc;


colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};

iLQR_Box_Time = cell(iScal,iRepts);
ExtMod_Box_Time = cell(iScal,iRepts);

iLQR_Box_Vbar = cell(iScal,iRepts);
ExtMod_Box_Vbar = cell(iScal,iRepts);


ScalerG = {};

ExtMod_normU = {}; iLQR_normU = {};
Exp_normU = {}; Tens_normU = {};
for idx = 1:iScal
    for i=1:iRepts
        try 
            iLQR_Box_Vbar{idx,i} = iLQR_Box{idx,i}.VnoMu{:}(end);     
            iLQR_Box_Time{idx,i} = iLQR_Box{idx,i}.Time;
        catch
            iLQR_Box_Vbar{idx,i} = nan; 
            iLQR_Box_Time{idx,i} = nan; 
        end
        try
            ExtMod_Box_Vbar{idx,i} =  ExtMod_Box{idx,i}.VnoMu{:}(end); 
            ExtMod_Box_Time{idx,i} =  ExtMod_Box{idx,i}.Time; 
        catch
            ExtMod_Box_Time{idx,i} = nan;
            ExtMod_Box_Vbar{idx,i} = nan; 
        end
            
        ScalerG{idx,i} = cell2mat(cWork{idx,i}(2));
        
%         Exp_Box_Time{idx,i}  = [Exp_Box{idx,i}.Time];     
%         Tens_Box_Time{idx,i} = [Tens_Box{idx,i}.Time];

%         iLQR_normU{idx,i} = norm(iLQR_Box{idx,i}.ubar{:});
%         ExtMod_normU{idx,i} = norm(ExtMod_Box{idx,i}.ubar{:});
%         Exp_normU{idx,i} = norm(Exp_Box{idx,i}.ubar{:});
%         Tens_normU{idx,i} = norm(Tens_Box{idx,i}.ubar{:});

    end
end


iLQR_Box_Vbar = cell2mat(iLQR_Box_Vbar); 
ExtMod_Box_Vbar = cell2mat(ExtMod_Box_Vbar);

iLQR_Box_Time = cell2mat(iLQR_Box_Time); 
ExtMod_Box_Time = cell2mat(ExtMod_Box_Time);

ScalerG = cell2mat(ScalerG);

%%
figure; hold on;
for idx = 1:iRepts
    h{1}=plot(ScalerG(:,idx),iLQR_Box_Vbar(:,idx),'ro','DisplayName','iLQR');
    h{2}=plot(ScalerG(:,idx),ExtMod_Box_Vbar(:,idx),'bo','DisplayName','DDP');
end
legend([h{:}])
xlabel('Noise Multiplicative','Interpreter','latex')
ylabel('End Cost, 20 Iterations','Interpreter','latex')
grid on; grid minor; 

%%
figure; 
hh = {names{1:size(Boxed_Time,1)}};
boxplot([iLQR_Box_Vbar(:) ExtMod_Box_Vbar(:)],hh,'Whisker',100);

DIST = 'normal';
figure; lw = 4;  nbins =10; light =0.2;

HH=histfit(iLQR_Box_Vbar(:),nbins,DIST);hold on; 
HH(2).DisplayName = names{1}; HH(2).LineWidth = lw; 

FD=histfit(ExtMod_Box_Vbar(:),nbins,DIST);hold on; 
FD(2).DisplayName = names{2}; FD(2).LineWidth = lw; 
legend

HH(1).FaceColor = colors{1}; HH(1).FaceAlpha = light; 
HH(2).Color = colors{1};
FD(1).FaceColor = colors{2}; FD(1).FaceAlpha = light;
FD(2).Color = colors{2};

Lgd =legend([HH(2),FD(2)]);
Lgd.Location = 'best';Lgd.FontSize  = 25;
set(gca,'FontSize',20);
% set(gca, 'XScale', 'log');
grid on; grid minor;

%%
figure; hold on;
for idx = 1:iRepts
    h{1}=plot(ScalerG(:,idx),iLQR_Box_Time(:,idx),'ro','DisplayName','iLQR');
    h{2}=plot(ScalerG(:,idx),ExtMod_Box_Time(:,idx),'bo','DisplayName','DDP');
end
legend([h{:}])
xlabel('Noise Multiplicative','Interpreter','latex')
ylabel('Time, 20 Iterations','Interpreter','latex')




%%
Boxed_Time = [mean(iLQR_Box_Time,'omitnan');mean(ExtMod_Box_Time,'omitnan')];
%     ...
%     Exp_Box_Time;Tens_Box_Time]';
figure; 
hh = {names{1:size(Boxed_Time,1)}};
boxplot(Boxed_Time',hh); 
grid on; grid minor
ylabel('Time(s), 20 Iterations','Interpreter','latex');
bb=gca;
bb.XAxis.TickLabelInterpreter='latex'; 
bb.XTickLabelRotation=15;
bb.YScale = 'log';






