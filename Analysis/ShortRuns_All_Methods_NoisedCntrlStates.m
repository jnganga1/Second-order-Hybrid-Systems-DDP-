%% Add noise to ubar - Run iLQR/DDP for 20 Iterations 
close all; clear; clc; 
% Add noise to ubar -- some weiner noise.
% Run iLQR/DDP for 50 Iterations and get time - makes figure A 
Itgr = 1;
OtherTests.Run = true; 
OtherTests.randomInit = false;
OtherTests.maxItr = 100;
OtherTests.dt = 1e-3;
OtherTests.BeforeNoise = false;


OtherTests.varyU  = true; %care about varying u in this round

iLQR_Box = {}; ExtMod_Box ={};
Exp_Box = {}; Tens_Box= {};

iNumber = 40; 

cWork = {};

%

for i = 1:iNumber    
%     N = length(ubar);
%
    %Control Noise
    N = 489; dt = 1e-3; 
    dW = sqrt(dt)*randn(4,N);   
    W = cumsum(dW) *  120;
    OtherTests.noise = W; 
    OtherTests.scalarNse =5;


    %State Noise
    StateScalar = 10;
    N = 489; dt = 1e-3; 
    dW = sqrt(dt)*randn(14,N);   
    W = cumsum(dW)*StateScalar(i); 
    W(8:end,:) = 0 *W(8:end,:);
    OtherTests.StateNoise = W;
    OtherTests.StateScalar = StateScalar(i);



    iLQR = 1; Method = 'none'; Reg=  1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method;
    OtherTests.RegMethod = 1;

    iLQR_Box{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR

    iLQR = 0; Method = 2; Reg = 1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 

    ExtMod_Box{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod

    iLQR = 0; Method = 1; Reg = 1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 

%     Exp_Box{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Exp
    iLQR = 0; Method = 3; Reg = 1;

    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 
%     Tens_Box{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Tens
end 

% figure; 
% plot(iLQR_Box{1}.Vbar{:},'DisplayName','iLQR'); hold on; 
% plot(ExtMod_Box{1}.Vbar{:},'DisplayName','ExtMod');
% plot(Exp_Box{1}.Vbar{:},'DisplayName','Exp');
% plot(Tens_Box{1}.Vbar{:},'DisplayName','Tens');
% legend

save('MatData/VaryU_Noise_ShortRuns_Not_Tens_fixed_ReModded2.mat');

%%
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)
close all; clear; clc;
% load('crcMatData/VaryU_Noise_MultiRuns_Again.mat');
% load('MatData/VaryU_Noise_ShortRuns_WithReg.mat');
% % load('MatData/VaryU_Noise_MultiRuns_Again.mat');
load('crcMatData_New/VaryU_Noise_ShortRuns_WithReg.mat');

%%
load('Data_BoundingNewNew/VaryU_Noise_ShortRuns_Not_Tens_fixed.mat');

%%

colors = {'r','b','k','c'};
names = {'iLQR', 'E-ModRNEA DDP','Explicit DDP','Tensor DDP'};

iLQR_Box_Time = cell(iScal,iRepts);
ExtMod_Box_Time = cell(iScal,iRepts);
Exp_Box_Time = cell(iScal,iRepts);
Tens_Box_Time = cell(iScal,iRepts);

iLQR_Box_Vbar = cell(iScal,iRepts);
ExtMod_Box_Vbar = cell(iScal,iRepts);
Exp_Box_Vbar = cell(iScal,iRepts);
Tens_Box_Vbar = cell(iScal,iRepts);

ScalerG = {};

iLQR_normU   = cell(iScal,iRepts);
ExtMod_normU = cell(iScal,iRepts);
Exp_normU    = cell(iScal,iRepts);
Tens_normU   = cell(iScal,iRepts);

for idx = 1:iScal
    for i=1:iRepts
        try 
            iLQR_Box_Vbar{idx,i} = iLQR_Box{idx,i}.VnoMu{:}(end);     
            iLQR_Box_Time{idx,i} = iLQR_Box{idx,i}.Time;
            iLQR_normU{idx,i}    = norm(iLQR_Box{idx,i}.ubar{:});
        catch
            iLQR_Box_Vbar{idx,i} = nan; 
            iLQR_Box_Time{idx,i} = nan; 
            iLQR_normU{idx,i} = nan;
        end
        try
            ExtMod_Box_Vbar{idx,i} =  ExtMod_Box{idx,i}.VnoMu{:}(end); 
            ExtMod_Box_Time{idx,i} =  ExtMod_Box{idx,i}.Time; 
            ExtMod_normU{idx,i}    = norm(ExtMod_Box{idx,i}.ubar{:});
        catch
            ExtMod_Box_Time{idx,i} = nan;
            ExtMod_Box_Vbar{idx,i} = nan; 
            ExtMod_normU{idx,i}  = nan; 
        end
        try
            Exp_Box_Vbar{idx,i} =  Exp_Box{idx,i}.VnoMu{:}(end); 
            Exp_Box_Time{idx,i} =  Exp_Box{idx,i}.Time; 
            Exp_normU{idx,i}    = norm(Exp_Box{idx,i}.ubar{:});
        catch
            Exp_Box_Time{idx,i} = nan;
            Exp_Box_Vbar{idx,i} = nan; 
            Exp_normU{idx,i} = nan; 
        end
        try
            Tens_Box_Vbar{idx,i} =  Tens_Box{idx,i}.VnoMu{:}(end); 
            Tens_Box_Time{idx,i} =  Tens_Box{idx,i}.Time; 
            Tens_normU{idx,i}    = norm(Tens_Box{idx,i}.ubar{:});
        catch
            Tens_Box_Time{idx,i} = nan;
            Tens_Box_Vbar{idx,i} = nan; 
            Tens_normU{idx,i} = nan; 
        end
        
        try            
        ScalerG{idx,i} = cell2mat(cWork{idx,i}(2));
        catch
            1==1; 
        end
    end
end


iLQR_Box_Vbar = cell2mat(iLQR_Box_Vbar); 
ExtMod_Box_Vbar = cell2mat(ExtMod_Box_Vbar);
Exp_Box_Vbar = cell2mat(Exp_Box_Vbar);
Tens_Box_Vbar = cell2mat(Tens_Box_Vbar);

iLQR_Box_Time = cell2mat(iLQR_Box_Time); 
ExtMod_Box_Time = cell2mat(ExtMod_Box_Time);
Exp_Box_Time = cell2mat(Exp_Box_Time);
Tens_Box_Time = cell2mat(Tens_Box_Time);


iLQR_normU   = cell2mat(iLQR_normU);
ExtMod_normU = cell2mat(ExtMod_normU);
Exp_normU    = cell2mat(Exp_normU);
Tens_normU   = cell2mat(Tens_normU);

ScalerG = cell2mat(ScalerG);

%%  Plot EndCost 
figure; hold on;
for idx = 1:iRepts
    h{1}=plot(ScalerG(:,idx),iLQR_Box_Vbar(:,idx),'ro','DisplayName',names{1},'Color',colors{1});
    h{2}=plot(ScalerG(:,idx),ExtMod_Box_Vbar(:,idx),'o','DisplayName',names{2},'Color',colors{2});
    h{3}=plot(ScalerG(:,idx),Exp_Box_Vbar(:,idx),'o','DisplayName',names{3},'Color',colors{3});
    h{4}=plot(ScalerG(:,idx),Tens_Box_Vbar(:,idx),'o','DisplayName',names{4},'Color',colors{4});


end
legend([h{:}])
xlabel('Noise Multiplicative','Interpreter','latex')
ylabel('End Cost, 20 Iterations','Interpreter','latex')
grid on; grid minor; 

%%  Boxplot EndCost 
figure; 
bp_vbar = boxplot([iLQR_Box_Vbar(:) ExtMod_Box_Vbar(:) ...
    Exp_Box_Vbar(:) Tens_Box_Vbar(:)],names,'Whisker',50);

ylabel('End Cost, 20 Iterations','Interpreter','latex')
bb=gca;
bb.XTickLabelRotation = 20; 
bb.TickLabelInterpreter = 'latex';
grid on; grid minor; 
%% Boxplot Time 
figure; 
boxplot([iLQR_Box_Time(:) ExtMod_Box_Time(:) ...
    Exp_Box_Time(:) Tens_Box_Time(:)],names,'Whisker',50);


gg=vpa(mean([iLQR_Box_Time(:) ExtMod_Box_Time(:) ...
    Exp_Box_Time(:) Tens_Box_Time(:)]),10)

ylabel('Time, 20 Iterations','Interpreter','latex')
bb=gca;
bb.YScale ='log';
bb.XTickLabelRotation = 20; 
bb.TickLabelInterpreter = 'latex';
grid on; grid minor; 


%% End Cost histfit -> This is pure binning technique, not probablity func
DIST = 'normal';
figure; lw = 4;  nbins =50; light =0.2;

iLQR_fit=histfit(iLQR_Box_Vbar(:),nbins,DIST);hold on; 
iLQR_fit(2).DisplayName = names{1}; iLQR_fit(2).LineWidth = lw; 

ExtMod_fit=histfit(ExtMod_Box_Vbar(:),nbins,DIST);hold on; 
ExtMod_fit(2).DisplayName = names{2}; ExtMod_fit(2).LineWidth = lw; 

Exp_fit=histfit(ExtMod_Box_Vbar(:),nbins,DIST);hold on; 
Exp_fit(2).DisplayName = names{3}; Exp_fit(2).LineWidth = lw; 

Tens_fit=histfit(ExtMod_Box_Vbar(:),nbins,DIST);hold on; 
Tens_fit(2).DisplayName = names{4}; Tens_fit(2).LineWidth = lw; 

%styling
legend

iLQR_fit(1).FaceColor = colors{1}; iLQR_fit(1).FaceAlpha = light; 
iLQR_fit(2).Color = colors{1};

ExtMod_fit(1).FaceColor = colors{2}; ExtMod_fit(1).FaceAlpha = light;
ExtMod_fit(2).Color = colors{2};

Exp_fit(1).FaceColor = colors{3}; Exp_fit(1).FaceAlpha = light;
Exp_fit(2).Color = colors{3};

Tens_fit(1).FaceColor = colors{4}; Tens_fit(1).FaceAlpha = light;
Tens_fit(2).Color = colors{4};

Lgd =legend([iLQR_fit(2),ExtMod_fit(2),Exp_fit(2),Tens_fit(2)]);
Lgd.Location = 'best';Lgd.FontSize  = 15;
set(gca,'FontSize',15);
% set(gca, 'XScale', 'log');
grid on; grid minor;
xlabel('End Cost','Interpreter','latex');
ylabel('Count - Normal Dist','Interpreter','latex');

%% End Cost bar/pdf -> both are pdfs

hh= {};
figure;
iLQR_y = iLQR_fit(1).YData;
iLQR_ybar = bar(iLQR_fit(1).XData, iLQR_y/sum(iLQR_y),colors{1});            % Probability
hold on
hh{1} = plot(iLQR_fit(2).XData, iLQR_fit(2).YData/sum(iLQR_y), colors{1},'DisplayName',names{1});


Exp_y = Exp_fit(1).YData;
Exp_ybar = bar(Exp_fit(1).XData, Exp_y/sum(Exp_y),colors{2});            % Probability
hh{2} = plot(Exp_fit(2).XData, Exp_fit(2).YData/sum(Exp_y), colors{2}, 'DisplayName',names{2});


ExtMod_y = ExtMod_fit(1).YData;
ExtMod_ybar = bar(ExtMod_fit(1).XData, ExtMod_y/sum(ExtMod_y),colors{3});            % Probability
hh{3} = plot(ExtMod_fit(2).XData, ExtMod_fit(2).YData/sum(ExtMod_y), colors{3},'DisplayName',names{3});

Tens_y = Tens_fit(1).YData;
Tens_ybar = bar(Tens_fit(1).XData, Tens_y/sum(Tens_y),colors{4});            % Probability
hh{4} = plot(Tens_fit(2).XData, Tens_fit(2).YData/sum(Tens_y), colors{4},'DisplayName',names{4});


iLQR_ybar.FaceAlpha = light;
ExtMod_ybar.FaceAlpha = light;
Exp_ybar.FaceAlpha = light;
Tens_ybar.FaceAlpha = light;

Lgd =legend([hh{:}]);
Lgd.Location = 'best';Lgd.FontSize  = 15;
set(gca,'FontSize',15);
% set(gca, 'XScale', 'log');
grid on; grid minor;
xlabel('End Cost','Interpreter','latex');
ylabel('Probability Density - Normal Dist','Interpreter','latex');
%% End Cost histogram bin count, pdf 
hh= {};
figure;
iLQR_y = iLQR_fit(1).YData;
% iLQR_ybar = bar(iLQR_fit(1).XData, iLQR_y/sum(iLQR_y),colors{1}); %Probability
yyaxis left
iLQR_histo = histogram(iLQR_Box_Vbar(:),'Normalization','probability','NumBins',nbins); hold on;
hold on
yyaxis right

hh{1} = plot(iLQR_fit(2).XData, iLQR_fit(2).YData/sum(iLQR_y), colors{1},'DisplayName',names{1});


%%
const = {'Probability Density - '}; 
str = {'Count','percentage','Normal Dist'};
normDist = {'count','probability', 'pdf'};

iCare = 1; 
nbins = 50;

m = 10^8;

hh = {}; 
figure; 
h = iLQR_Box_Vbar; iClr = 1;
yyaxis left
iLQR_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;
xlabel('End Cost','Interpreter','latex'); 
ylabel(strcat(const,str{iCare}),'Interpreter','latex');

yyaxis right
[L,D]= ksdensity(h(:),'Function','pdf');
hh{iClr} = plot(D,L*m,colors{iClr},'DisplayName',names{iClr});  
ylabel('Probablity Density','Interpreter','latex');
bb=gca; 


%
% figure;
h = Exp_Box_Vbar; iClr = 2;
yyaxis left
Exp_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;

yyaxis right
[L,D]= ksdensity(h(:),'Function','pdf');
hh{iClr} = plot(D,L*m,colors{iClr},'DisplayName',names{iClr});  
%
%
yyaxis left 
h = ExtMod_Box_Vbar; iClr = 3;
ExtMod_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;

yyaxis right 
[L,D]= ksdensity(h(:),'Function','pdf');
hh{iClr} = plot(D,L*m,colors{iClr},'DisplayName',names{iClr});  


yyaxis left
h = Tens_Box_Vbar; iClr = 4;
Tens_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;

yyaxis right
[L,D]= ksdensity(h(:),'Function','pdf');
hh{iClr} = plot(D,L*m,colors{iClr},'DisplayName',names{iClr});  

for i=1:4 
    hh{i}.LineStyle = '-';
    hh{i}.LineWidth = 2;
end

iLQR_histo.FaceColor = colors{1};
Exp_histo.FaceColor = colors{2};
ExtMod_histo.FaceColor = colors{3};
Tens_histo.FaceColor = colors{4};

iLQR_histo.FaceAlpha = light; 
Exp_histo.FaceAlpha = light; 
ExtMod_histo.FaceAlpha = light; 
Tens_histo.FaceAlpha = light; 


Lgd =legend([hh{:}]);
Lgd.Location = 'best';Lgd.FontSize  = 15;
set(gca,'FontSize',15);
% set(gca, 'XScale', 'log');
grid on; grid minor;
xlabel('End Cost','Interpreter','latex');
ylabel('Probability Density - Normal Dist','Interpreter','latex');

%%

%%
%{
rng('default'); 
x = 100*[randn(30,1); 5+randn(30,1)];
[f,xi,bw] = ksdensity(x); 
figure
plot(xi,f*10^2);
xlabel('xi')
ylabel('f')
hold on
%}
%






