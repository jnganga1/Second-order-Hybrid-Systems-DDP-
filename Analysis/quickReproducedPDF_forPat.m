set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',22)



DDP = h.Children(1).Data;
iLQR = h.Children(2).Data;

iLQR(iLQR > 2.5*10^7) =[]; 
DDP(DDP > 2.5*10^7) =[]; 


figure;
yyaxis left 

g{1} = histogram(iLQR);
hold on
g{2} = histogram(DDP);

clr = {'r','b'};
names = {'DDP','iLQR'};
for i=1:2 
    g{i}.BinWidth =263000; 
    g{i}.DisplayName = names{i}; 
    g{i}.FaceAlpha =0.4000; 
    g{i}.FaceColor = clr{i};    
end
xlabel('End Cost to Optimality','Interpreter','latex');
ylabel('Count','Interpreter','latex');


pd_fits = {}; 
pdfs_all = {};
linePlots = {}; 
pdf_val = 'Normal'; 

x_values = -2e7:1e3:3e7;
m=1;

yyaxis right 
All_vbars = {DDP,iLQR};
for i = 1:length(All_vbars) 
    pd_fits{i} = fitdist(All_vbars{i},pdf_val);
    pdfs_all{i} = pdf( pd_fits{i},x_values);
    linePlots{i}=plot(x_values,pdfs_all{i}*m,'LineWidth',2,'DisplayName',names{i},...
        'Color',clr{i},'LineStyle','-'); hold on
    fprintf('%s area under curve: %.2f\n', names{i},trapz(x_values,pdfs_all{i}));
end
m  = get(gca,'Children');
set(gca,'Children',[flipud(m)])
legend([linePlots{:}]);
grid on; grid minor;



% legend