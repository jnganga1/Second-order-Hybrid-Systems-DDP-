

% 
ilqr=[]; tens=[];exp=[];extmod=[];
numRepts = 2;
for irepts = 1:numRepts
    iLQR = 1; Method = 'none';
    [Time,ilqr_Vstore,ilqr_norm_hbar]  = DDP_2DQuadruped(iLQR,Method,0);
    ilqr(end+1) = Time; 
    
    %1 is Exp, 2 is ExtMod, 3 is Tensor
    iLQR = 0; Method = 1;
    [Time,exp_Vstore,exp_norm_hbar] =DDP_2DQuadruped(iLQR,Method);
    exp(end+1) = Time;
    
    iLQR = 0; Method = 2; regularizationMethod = 1;
    [Time,extmod_Vstore,extmod_norm_hbar] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);
    extmod(end+1) = Time;
    
    iLQR = 0; Method = 3;
   [Time,tens_Vstore,tens_norm_hbar]  =DDP_2DQuadruped(iLQR,Method);  
   tens(end+1) = Time;
end
1==1;




set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

figure;
labels = {'iLQR', 'Explicit Method', 'Extended RNEA', 'Tensor'};
colored = {'b','m','g','r'};
all = [ilqr' exp' extmod' tens'];
Tme = boxplot(all,'labels',labels); %'whisker',2, 'widths',0.4,'symbol','');
set(gca, 'YScale', 'log'); 
set(gca,'FontSize',30);
set(Tme,'LineWidth', 2);
L = gca;
L.XAxis.TickLabelInterpreter = 'latex';
grid on; grid minor; 
xtickangle(15);
title('DDP Evaluation (3 repetitions)')
1==1;

figure; lw =2;
plot(ilqr_norm_hbar,'LineWidth',lw,'Color',colored{1}); hold on
plot(exp_norm_hbar,'LineWidth',lw,'Color',colored{2})
plot(extmod_norm_hbar,'LineWidth',lw,'Color',colored{3})
plot(tens_norm_hbar,'LineWidth',lw,'Color',colored{4})
set(gca,'FontSize',30);
legend(labels)


figure; 
semilogy(ilqr_Vstore - ilqr_Vstore(end),'LineWidth',lw,'Color',colored{1}); hold on 
semilogy(exp_Vstore - exp_Vstore(end),'LineWidth',lw,'Color',colored{2});
semilogy(extmod_Vstore - extmod_Vstore(end),'LineWidth',lw,'Color',colored{3});
semilogy(tens_Vstore - tens_Vstore(end),'LineWidth',lw,'Color',colored{4});
set(gca,'FontSize',30);
legend(labels)
grid on; grid minor