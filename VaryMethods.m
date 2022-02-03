
%Using iLQR method - no regularization
iLQR = 1; Method = 'none';
[ilqr_Time, ilqr_Vstore,ilqr_xbar,ilqr_norm_hbar,...
    ilqr_ybar,ilqr_ubar,ilqr_VstoreNoMu]  = DDP_2DQuadruped(iLQR,Method,0);

% %using Exp Regularization Method
iLQR = 0; Method = 1; regularizationMethod = 1;
[Exp_Time,Exp_Vstore,Exp_Xbar,Exp_norm_hbar,...
    Exp_ybar,Exp_ubar,Exp_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);

% %using ExtMod Regularization Method
iLQR = 0; Method = 2; regularizationMethod = 1;
[ExtMod_Time,ExtMod_Vstore,ExtMod_Xbar,ExtMod_norm_hbar,...
    ExtMod_ybar,ExtMod_ubar,ExtMod_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);

% %using Traditional Regularization Method
iLQR = 0; Method = 3; regularizationMethod = 1;
[TradReg_Time,TradReg_Vstore,TradReg_norm_Xbar,TradReg_norm_hbar,...
    TradReg_ybar,TradReg_ubar,TradReg_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);

% load('TheMethods_for_Wed.mat')
%%
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

maxItr = 100;

colors ={'c','r','b','k'};
figure; subplot(2,1,1);
% yyaxis left
plot(ilqr_VstoreNoMu,colors{1},'DisplayName','iLQR'); 
hold on; 
plot(Exp_VstoreNoMu ,colors{2},'DisplayName','Explicit Method');
plot(ExtMod_VstoreNoMu,colors{3},'DisplayName','Extended Modified Method');
plot(TradReg_VstoreNoMu,colors{4},'DisplayName','Tensor Method');
set(gca,'Yscale','log')
% title('Cost (log)');
ylabel('Cost w/o constraint (log)'); 
xlabel('Iterations');

%
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


% yyaxis right
subplot(2,1,2);
colors ={'c','r','b','k'};
txt = ['iLQR, Time: = ',num2str(ilqr_Time)];
plot(ilqr_norm_hbar,colors{1},'DisplayName',txt);

hold on;
txt = ['Explicit Method,Time: = ',num2str(Exp_Time)];
plot(Exp_norm_hbar,colors{2},'DisplayName',txt)

txt = ['Extended Modified, Time: = ',num2str(ExtMod_Time)];
plot(ExtMod_norm_hbar,colors{3},'DisplayName',txt)

txt = ['Tensor Method, Time: =' num2str(TradReg_Time)];
plot(TradReg_norm_hbar,colors{4},'DisplayName',txt)
set(gca,'Yscale','log'); 
hold all

ylabel('Constraint Violation'); 
xlabel('Iterations');

x = [1 1]*maxItr; 
h = gca; hold on 
for i= 1:length(ilqr_Vstore)/100
    l=plot(i*x,h.YLim,'k','LineWidth',2); 
end
hold off
h = legend;
h.Location = 'best';
h.String(5:end) = '';
% title('Impact Event contact violation');



%%
set(groot,'defaultLineLineWidth',1)
figure; 
forces = {'Torque,hip,front','Torque,knee,front',...
    'Torque,hip,back','Torque,hip,back'};
time = linspace(0,2.93,293);
for i=1:4
subplot(2,2,i);
plot(ilqr_ubar(i,:),colors{1},'DisplayName','iLQR'); hold on;
plot(Exp_ubar(i,:),colors{2},'DisplayName','Explicit Method');
plot(ExtMod_ubar(i,:),colors{3},'DisplayName','Extended Modified Method');
plot(TradReg_ubar(i,:),colors{4},'DisplayName','Tensor Method');
ylabel(forces{i});
% xlabel('Time')
end
set(gca,'LineWidth',1)
h =legend;
h.Location = 'best';
%%
i=2;
figure;
plot(ilqr_ubar(i,:),colors{1},'DisplayName','iLQR'); hold on;
plot(TradReg_ubar(i,:),colors{4},'DisplayName','Tensor Method');

