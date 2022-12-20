%assumes right data loaded 
%%
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',22)
close all; clearvars; clc;
%%

colors = {'r','b','k','c'};
names = {'iLQR', 'mRNEAc DDP','Tensor DDP'};

iLQR_Box_Time = cell(iNumber);
ExtMod_Box_Time = cell(iNumber);
Exp_Box_Time = cell(iNumber);
Tens_Box_Time = cell(iNumber);

iLQR_Box_Vbar = cell(iNumber);
ExtMod_Box_Vbar = cell(iNumber);
Exp_Box_Vbar = cell(iNumber);
Tens_Box_Vbar = cell(iNumber);

ScalerG = {};

iLQR_normU   = cell(iNumber);
ExtMod_normU = cell(iNumber);
Exp_normU    = cell(iNumber);
Tens_normU   = cell(iNumber);

for idx = 1:iNumber
        try 
            iLQR_Box_Vbar{idx} = iLQR_Box{idx}.VnoMu{:}(end);     
            iLQR_Box_Time{idx} = iLQR_Box{idx}.Time;
            iLQR_normU{idx}    = norm(iLQR_Box{idx}.ubar{:});
        catch
            iLQR_Box_Vbar{idx} = nan; 
            iLQR_Box_Time{idx} = nan; 
            iLQR_normU{idx} = nan;
        end
        try
            ExtMod_Box_Vbar{idx} =  ExtMod_Box{idx}.VnoMu{:}(end); 
            ExtMod_Box_Time{idx} =  ExtMod_Box{idx}.Time; 
            ExtMod_normU{idx}    = norm(ExtMod_Box{idx}.ubar{:});
        catch
            ExtMod_Box_Time{idx} = nan;
            ExtMod_Box_Vbar{idx} = nan; 
            ExtMod_normU{idx}  = nan; 
        end
%         try
%             Exp_Box_Vbar{idx} =  Exp_Box{idx}.VnoMu{:}(end); 
%             Exp_Box_Time{idx} =  Exp_Box{idx}.Time; 
%             Exp_normU{idx}    = norm(Exp_Box{idx}.ubar{:});
%         catch
%             Exp_Box_Time{idx} = nan;
%             Exp_Box_Vbar{idx} = nan; 
%             Exp_normU{idx} = nan; 
%         end
        try
            Tens_Box_Vbar{idx} =  Tens_Box{idx}.VnoMu{:}(end); 
            Tens_Box_Time{idx} =  Tens_Box{idx}.Time; 
            Tens_normU{idx}    = norm(Tens_Box{idx}.ubar{:});
        catch
            Tens_Box_Time{idx} = nan;
            Tens_Box_Vbar{idx} = nan; 
            Tens_normU{idx} = nan; 
        end
        
%         try            
        ScalerG{idx} = iLQR_Box{idx}.OtherTests.scalarNse;
%         catch
%             1==1; 
%         end
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
h={};
figure; hold on;
% for idx = 1:iNumber
    h{1}=plot(iLQR_Box_Vbar,'ro','DisplayName',names{1},'Color',colors{1});
    h{2}=plot(ExtMod_Box_Vbar,'o','DisplayName',names{2},'Color',colors{2});
%     h{3}=plot(Exp_Box_Vbar,'o','DisplayName',names{3},'Color',colors{3});
    h{4}=plot(Tens_Box_Vbar,'o','DisplayName',names{3},'Color',colors{4});


% end
legend([h{:}])
xlabel('Noise Multiplicative','Interpreter','latex')
ylabel('End Cost, 50 Iterations','Interpreter','latex')
grid on; grid minor; 
legend
%%  Boxplot EndCost 
figure; 

bp_vbar = boxplot([iLQR_Box_Vbar(:) ExtMod_Box_Vbar(:)])% ...
%    Tens_Box_Vbar(:)],names,'Whisker',50);

% bp_vbar = boxplot([iLQR_Box_Vbar(:) ExtMod_Box_Vbar(:)]',names,'Whisker',50);

ylabel('End Cost, 50 Iterations','Interpreter','latex')
bb=gca;
bb.XTickLabelRotation = 20; 
bb.TickLabelInterpreter = 'latex';
grid on; grid minor; 
%% Boxplot Time 
figure; 
divy=1;
boxplot([iLQR_Box_Time(:) ExtMod_Box_Time(:) ...
    Exp_Box_Time(:) Tens_Box_Time(:)]./divy,names,'symbol', '');


colored = {'r','m','b'};
h = findobj(gca,'Tag','Box');
g = findobj(gca,'Tag','Median');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colored{j},'FaceAlpha',1);
    g(j).Color = 'w';
    g(j).LineWidth = 1;
end
set(gca,'children',flipud(get(gca,'children')))
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');


% boxplot([iLQR_Box_Time(:) ExtMod_Box_Time(:)]./divy,names,'Whisker',1000);

gg=vpa(mean([iLQR_Box_Time(:) ExtMod_Box_Time(:) ...
    Tens_Box_Time(:)]/divy),10)

for i = 1:2 
    gg(i+1)/gg(1)
end

ylabel('Time (s), 50 Iterations','Interpreter','latex')
bb=gca;
bb.YTick = [10^2 10^3];
bb.YScale ='log';
bb.XTickLabelRotation = 20; 
bb.TickLabelInterpreter = 'latex';
grid on; grid minor; 

%%
figure; hold on 
% hold on;
% rng default; % For reproducibility
r = normrnd(10,1,100,1);

pd = fitdist(r,'Normal');

r=histfit(r)
 trapz(r(2).XData,r(2).YData)
 
pd_ilqr=fitdist(iLQR_Box_Vbar(:),'LogNormal');
pd_ddp=fitdist(ExtMod_Box_Vbar(:),'LogNormal'); 

m=10^5; 


figure;
yyaxis right

x_values = -1e7:1e3:3e7;
y_ilqr = pdf(pd_ilqr,x_values);
y_ddp = pdf(pd_ddp,x_values);


plot(x_values,y_ilqr*m,'LineWidth',2,'DisplayName','iLQR'); hold on
plot(x_values,y_ddp*m,'LineWidth',2,'DisplayName','DDP'); 
legend
% h=gca;h=h.Children;set(gca,'Children',[flipud(h)])

%%
figure; 
All_vbars = { iLQR_Box_Vbar,ExtMod_Box_Vbar};
%     Exp_Box_Vbar,Tens_Box_Vbar};

data = load('MatData/NowOptData');
Vopt = data.ExtMod_Nsd_nowOpt.Vbar{1}(end);

pd_fits = {}; 
pdfs_all = {};
linePlots = {}; 
pdf_val = 'Normal'; 

x_values = -2e7:1e3:3e7;
m=1;
% m = 10^5; 

yyaxis right 
for i = 1:length(All_vbars) 
    pd_fits{i} = fitdist(All_vbars{i} - Vopt,pdf_val);
    pdfs_all{i} = pdf( pd_fits{i},x_values);
    linePlots{i}=plot(x_values,pdfs_all{i}*m,'LineWidth',2,'DisplayName',names{i},...
        'Color',colors{i},'LineStyle','-'); hold on
    fprintf('%s area under curve: %.2f\n', names{i},trapz(x_values,pdfs_all{i}));
end
h = get(gca,'Children');
set(gca,'Children',[flipud(h)])
legend([linePlots{:}]);
grid on; grid minor;


yyaxis left
nbins = 70;

hist = {};
for idx = 1:length(All_vbars)
   hist{idx} = histogram(All_vbars{idx}-Vopt,nbins,'DisplayName',names{idx}); hold on; 
   hist{idx}.FaceColor =  colors{idx};
   hist{idx}.FaceAlpha = 0.4;
end
xlabel('End Cost to Optimality','Interpreter','latex');
ylabel('Count','Interpreter','latex');

legend([linePlots{:}]);


%% End Cost histfit -> This is pure binning technique, not probablity func
DIST = 'gamma';
figure; lw = 4;  nbins =100; light =0.2;

iLQR_fit=histfit(iLQR_Box_Vbar(:),nbins,DIST);hold on; 
iLQR_fit(2).DisplayName = names{1}; iLQR_fit(2).LineWidth = lw; 

ExtMod_fit=histfit(ExtMod_Box_Vbar(:),nbins,DIST);hold on; 
ExtMod_fit(2).DisplayName = names{2}; ExtMod_fit(2).LineWidth = lw; 

% Exp_fit=histfit(ExtMod_Box_Vbar(:),nbins,DIST);hold on; 
% Exp_fit(2).DisplayName = names{3}; Exp_fit(2).LineWidth = lw; 

Tens_fit=histfit(ExtMod_Box_Vbar(:),nbins,DIST);hold on; 
Tens_fit(2).DisplayName = names{3}; Tens_fit(2).LineWidth = lw; 

%

%styling
legend

iLQR_fit(1).FaceColor = colors{1}; iLQR_fit(1).FaceAlpha = light; 
iLQR_fit(2).Color = colors{1};

ExtMod_fit(1).FaceColor = colors{2}; ExtMod_fit(1).FaceAlpha = light;
ExtMod_fit(2).Color = colors{2};

% Exp_fit(1).FaceColor = colors{3}; Exp_fit(1).FaceAlpha = light;
% Exp_fit(2).Color = colors{3};

Tens_fit(1).FaceColor = colors{3}; Tens_fit(1).FaceAlpha = light;
Tens_fit(2).Color = colors{3};

Lgd =legend([iLQR_fit(2),ExtMod_fit(2),Tens_fit(2)]);
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


% Exp_y = Exp_fit(1).YData;
% Exp_ybar = bar(Exp_fit(1).XData, Exp_y/sum(Exp_y),colors{2});            % Probability
% hh{2} = plot(Exp_fit(2).XData, Exp_fit(2).YData/sum(Exp_y), colors{2}, 'DisplayName',names{2});
% 

ExtMod_y = ExtMod_fit(1).YData;
ExtMod_ybar = bar(ExtMod_fit(1).XData, ExtMod_y/sum(ExtMod_y),colors{3});            % Probability
hh{3} = plot(ExtMod_fit(2).XData, ExtMod_fit(2).YData/sum(ExtMod_y), colors{3},'DisplayName',names{2});

Tens_y = Tens_fit(1).YData;
Tens_ybar = bar(Tens_fit(1).XData, Tens_y/sum(Tens_y),colors{4});            % Probability
hh{4} = plot(Tens_fit(2).XData, Tens_fit(2).YData/sum(Tens_y), colors{4},'DisplayName',names{3});


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

iCare = 2; 
nbins = 100;

m = 1;

hh = {}; 
figure; 
h = iLQR_Box_Vbar; iClr = 1;
yyaxis left
iLQR_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;
xlabel('End Cost','Interpreter','latex'); 
ylabel(strcat(const,str{iCare}),'Interpreter','latex');

yyaxis right
[L,D]= ksdensity(h(:),'Function','pdf');
hh{iClr} = plot(D,colors{iClr},'DisplayName',names{iClr});  
ylabel('Probablity Density','Interpreter','latex');
bb=gca; 
fprintf('iLQR trapz: %.2f\n',trapz(L,D));

% %
% % figure;
% h = Exp_Box_Vbar; iClr = 2;
% yyaxis left
% Exp_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;
% 
% yyaxis right
% [L,D]= ksdensity(h(:),'Function','pdf');
% hh{iClr} = plot(D,L*m,colors{iClr},'DisplayName',names{iClr});
% fprintf('Exp trapz: %.2f\n',trapz(L,D));
%
%
yyaxis left 
h = ExtMod_Box_Vbar; iClr = 2;
ExtMod_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;

yyaxis right 
[L,D]= ksdensity(h(:),'Function','pdf');
hh{iClr} = plot(D,L*m,colors{iClr},'DisplayName',names{iClr});  

% %
% yyaxis left
% h = Tens_Box_Vbar; iClr = 3;
% Tens_histo = histogram(h(:),'Normalization',normDist{iCare},'NumBins',nbins); hold on;
% 
% yyaxis right
% [L,D]= ksdensity(h(:),'Function','pdf');
% hh{iClr} = plot(D,L*m,colors{iClr},'DisplayName',names{iClr});  

% for i=1:4 
%     hh{i}.LineStyle = '-';
%     hh{i}.LineWidth = 2;
% end

iLQR_histo.FaceColor = colors{1};
% Exp_histo.FaceColor = colors{2};
ExtMod_histo.FaceColor = colors{3};
% Tens_histo.FaceColor = colors{4};

iLQR_histo.FaceAlpha = light; 
% Exp_histo.FaceAlpha = light; 
ExtMod_histo.FaceAlpha = light; 
% Tens_histo.FaceAlpha = light; 

%
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






