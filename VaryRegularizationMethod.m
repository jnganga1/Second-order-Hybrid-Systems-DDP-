
%Using iLQR method - no regularization
iLQR = 1; Method = 'none';
[ilqr_Time, ilqr_Vstore,ilqr_xbar,ilqr_norm_hbar,...
    ilqr_ybar,ilqr_ubar,ilqr_VstoreNoMu,iqr_store]  = DDP_2DQuadruped(iLQR,Method,0);
%%

%using TensorPullback Regularization Method
% iLQR = 0; Method = 2; regularizationMethod = 3;
% [Pullback_Time,Pullback_Vstore,Pullback_xbar,Pullback_norm_hbar,...
%     Pullback_ybar,Pullback_ubar,Pullback_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);

%
% %using Traditional Regularization Method

iLQR = 0; Method = 2; regularizationMethod = 1;
[TradReg_Time,TradReg_Vstore,TradReg_norm_Xbar,TradReg_norm_hbar,...
    TradReg_ybar,TradReg_ubar,TradReg_VstoreNoMu,DDP_store] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);

 save('TorquesConstraint.mat'); 
1==1;
%%
%using point-by-point Regularization Method
iLQR = 0; Method = 2; regularizationMethod = 2;
[Point_Time,Point_Vstore,Point_xbar,Point_norm_hbar,...
    Point_ybar,Point_ubar,Point_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);
save('WarmStart_RegComparison.mat');
%%

%using TensorPullback Regularization Method
iLQR = 0; Method = 2; regularizationMethod = 3;
[Pullback_Time,Pullback_Vstore,Pullback_xbar,Pullback_norm_hbar,...
    Pullback_ybar,Pullback_ubar,Pullback_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);


% %using Traditional Regularization Method
iLQR = 0; Method = 1; regularizationMethod = 1;
[TradRegMethod1_Time,TradRegMethod1_Vstore,TradRegMethod1_norm_Xbar,TradRegMethod1_norm_hbar,...
    TradRegMethod1_ybar,TradRegMethod1_ubar,TradRegMethod1_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);

% %using Traditional Regularization Method
iLQR = 0; Method = 3; regularizationMethod = 1;
[TradRegMethod3_Time,TradRegMethod3_Vstore,TradRegMethod3_norm_Xbar,TradRegMethod3_norm_hbar,...
    TradRegMethod3_ybar,TradRegMethod3_ubar,TradRegMethod3_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);



save('PleaseWorkAgain_Yep.mat');
%%
save('DoesNewWork.mat')

%%
fname = strcat('ConstraintAdded_ilQR ','.gif');
params.filename = fname;
simulation_Jump(ilqr_xbar(1:params.q_size,:),robot_params,params,true);


%%
% %NewMethods
% iLQR = 0; Method = 1; regularizationMethod = 2;
% [Exp_Time,Exp_Vstore,Exp_norm_hbar,...
%     Exp_ybar,Exp_ubar,Exp_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);
% 
% iLQR = 0; Method = 3; regularizationMethod = 2;
% [Tens_Time,Tens_Vstore,Tens_norm_hbar,...
%     Tens_ybar,Tens_ubar,Tens_VstoreNoMu] =DDP_2DQuadruped(iLQR,Method,regularizationMethod);


% save('RegularizationData_Better3.mat')
%
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

count_care = 4; 
count_care = count_care + 1; 
maxItr = 150;

colors ={'c','r','b','k'};
figure; subplot(2,1,1);
% yyaxis left
plot(ilqr_VstoreNoMu,colors{1},'DisplayName','iLQR - no regularization'); 
hold on; 
plot(TradReg_VstoreNoMu ,colors{2},'DisplayName','Traditional Regularization');
% plot(Point_VstoreNoMu,colors{3},'DisplayName','Time-point Regularization');
% plot(Pullback_VstoreNoMu,colors{4},'DisplayName','Tensor Term Scaleback Regularization');
set(gca,'Yscale','log')
% title('Cost (log)');
ylabel('Cost w/o constraint (log)'); 
xlabel('Iterations'); grid on; grid minor;

hold all
x = [1 1]*maxItr; 
h = gca; hold on 
for i= 1:length(ilqr_Vstore)/maxItr
    l=plot(i*x,h.YLim,'k','LineWidth',2); 
end
hold off
h = legend;
h.Location = 'best';
h.String(count_care:end) = '';
%
% yyaxis right
subplot(2,1,2);
%%
figure;
colors ={'c','r','b','k'};
txt = ['iLQR - no regularization, Time: = ',num2str(ilqr_Time)];
plot(ilqr_norm_hbar,colors{1},'DisplayName',txt);
% 
hold on;
% txt = ['DDP - Time-point regularization, Time: = ',num2str(Point_Time)];
% plot(Point_norm_hbar,colors{2},'DisplayName',txt)

txt = ['DDP - Traditional regularization, Time: = ',num2str(TradReg_Time)];
plot(TradReg_norm_hbar,colors{3},'DisplayName',txt)

% txt = ['DDP - Tensor Term Pullback regularization, Time: =' num2str(Pullback_Time)];
% plot(Pullback_norm_hbar,colors{4},'DisplayName',txt)
set(gca,'Yscale','linear'); 
hold all
grid on; grid minor;
%
ylabel('Constraint Violation'); 
xlabel('Iterations');
%%
x = [1 1]*maxItr; 
h = gca; hold on 
for i= 1:length(ilqr_Vstore)/maxItr
    l=plot(i*x,h.YLim,'k','LineWidth',2); 
end
hold off
% h = legend;
% h.Location = 'best';
% h.String(count_care:end) = '';
% title('Impact Event contact violation');
%
set(groot,'defaultLineLineWidth',1)
figure; 
forces = {'Torque,hip,front','Torque,knee,front',...
    'Torque,hip,back','Torque,hip,back'};
%%
time = linspace(0,2.93,293);
for i=1:4
subplot(2,2,i);
plot(ilqr_ubar(i,:),colors{1},'DisplayName','iLQR - no regularization'); hold on;
plot(TradReg_ubar(i,:),colors{2},'DisplayName','Traditional Regularization');
% plot(Point_ubar(i,:),colors{3},'DisplayName','Time-point Regularization');
% plot(Pullback_ubar(i,:),colors{4},'DisplayName','Tensor Term Scaleback Regularization');
ylabel(forces{i});
grid on; grid minor; 
h=gca; h_len = length(h.Children(1).XData); 
xlim([0 h_len]);
% xlabel('Time')
end
set(gca,'LineWidth',1)
h =legend;
h.Location = 'best';
%
%
set(groot,'defaultLineLineWidth',4)
figure; 
forces = {'fx, front','fz,front'};

time = linspace(0,2.93,293);
%%
figure;
%Front
for i=1:2
    subplot(2,2,i); hold on;
    plot(ilqr_ybar(i,:,1),colors{1},'DisplayName','iLQR - no regularization');
    plot(TradReg_ybar(i,:,1),colors{2},'DisplayName','Traditional Regularization');
%     plot(Point_ybar(i,:,1),colors{3},'DisplayName','Time-point Regularization');
%     plot(Pullback_ybar(i,:,1),colors{4},'DisplayName','Tensor Term Pullback Regularization');
    ylabel(forces{i});
    grid on; grid minor; 
    h=gca; h_len = length(h.Children(1).XData); 
    xlim([0 h_len]);
    xlabel('Time')
end

%
forces = {'fx, back','fz,back'};
time = linspace(0,2.93,293);

for i=1:2
    subplot(2,2,2+i); hold on;
    plot(ilqr_ybar(i,:,2),colors{1},'DisplayName','iLQR - no regularization');
    plot(TradReg_ybar(i,:,2),colors{2},'DisplayName','Traditional Regularization');
    plot(Point_ybar(i,:,2),colors{3},'DisplayName','Time-point Regularization');
    plot(Pullback_ybar(i,:,2),colors{4},'DisplayName','Tensor Term Scaleback Regularization');
    ylabel(forces{i});
    grid on; grid minor; 
    h=gca; h_len = length(h.Children(1).XData); 
    xlim([0 h_len]);
    xlabel('Time')
end
set(gca,'LineWidth',.05)
h =legend;
h.Location = 'best';

set(gca,'LineWidth',1)
h =legend;
h.Location = 'best';

%
colors ={'c','r','b','k'};
figure; 
% yyaxis left
plot(ilqr_Vstore,colors{1},'DisplayName','iLQR - no regularization'); 
hold on; 
plot(TradReg_Vstore ,colors{2},'DisplayName','Traditional Regularization');
plot(Point_Vstore,colors{3},'DisplayName','Time-point Regularization');
plot(Pullback_Vstore,colors{4},'DisplayName','Tensor Term Scaleback Regularization');
set(gca,'Yscale','log')
% title('Cost (log)');
ylabel('Cost w/ constraint term (log)'); 
xlabel('Iterations');
%
hold all
x = [1 1]*maxItr; 
h = gca; hold on 
for i= 1:length(ilqr_Vstore)/maxItr 
    l=plot(i*x,h.YLim,'k','LineWidth',2); 
end
hold off
h = legend;
h.Location = 'best';
h.String(count_care:end) = '';

%%

set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

colors ={'c','r','b','k'};

k_AL = 5; 
maxItr_space = 0:maxItr:maxItr*k_AL; maxItr_space(1) =1;
figure;
for i=1:k_AL
    begin = maxItr_space(i)+1; last = maxItr_space(i+1); 
    semilogy(begin:last, ilqr_Vstore(begin:last) - ilqr_Vstore(last)...
        ,colors{1},'DisplayName','iLQR - no regularization'); hold on;
%     semilogy(begin:last,Pullback_Vstore(begin:last) - Pullback_Vstore(last)...
%         ,colors{4},'DisplayName','Tensor Term scaleback Regularization');
    semilogy(begin:last,TradReg_Vstore(begin:last)-TradReg_Vstore(last)...
        ,colors{2},'DisplayName','Traditional Regularization');
%     semilogy(begin:last, Point_Vstore(begin:last) - Point_Vstore(last)...
%         ,colors{3},'DisplayName','Time-point Regularization');
end
%
grid on; grid minor; 
g= gca;
h = legend;
h.Location = 'best';
h.String(count_care:end) ='';

hold all
x = [1 1]*maxItr; 
h = gca; hold on 
for i= 1:length(ilqr_Vstore)/maxItr 
    l=plot(i*x,h.YLim,'k','LineWidth',2); 
end
hold off
h = legend;
h.Location = 'best';
h.String(count_care:end) = '';

%%
colors ={'c','r','b','k'};
figure; 
% yyaxis left
semilogy(ilqr_Vstore-ilqr_Vstore(end),colors{1},'DisplayName','iLQR - no regularization'); 
hold on; 
semilogy(TradReg_Vstore-TradReg_Vstore(end),colors{2},'DisplayName','Traditional Regularization');
% semilogy(Point_Vstore-Point_Vstore(end),colors{3},'DisplayName','Time-point Regularization');
% semilogy(Pullback_Vstore-Pullback_Vstore(end),colors{4},'DisplayName','Tensor Term Scaleback Regularization');
% set(gca,'Yscale','log')
% title('Cost (log)');
ylabel('Cost w/ constraint term (log)'); 
xlabel('Iterations');

hold all
x = [1 1]*maxItr; 
h = gca; hold on 
for i= 1:length(ilqr_Vstore)/maxItr 
    l=plot(i*x,h.YLim,'k','LineWidth',2); 
end
hold off
h = legend;
h.Location = 'best';
h.String(count_care:end) = '';


%%  
figure;
k = 0; 
for i =  1:length(iqr_store.Vbar) 
    l = length(iqr_store.Vbar{i}); 
    semilogy(k:k+l-1,iqr_store.Vbar{i} - iqr_store.Vbar{i}(end),'c'); hold on; 
    k = k +l+1; 
end
hold on; k = 0; 
for i =  1:length(DDP_store.Vbar) 
    l = length(DDP_store.Vbar{i}); 
    semilogy(k:k+l-1,DDP_store.Vbar{i} - DDP_store.Vbar{i}(end),'r'); hold on; 
    k = k +l+1; 
end
grid on; grid minor; 
xlabel('Iterations'); 
ylabel('Value'); 
h.XLabel.Interpreter = 'latex';
h.YLabel.Interpreter= 'latex';

%%
figure;
k = 0; 
for i =  1:length(iqr_store.Vbar) 
    semilogy(k:k+l-1,iqr_store.Vbar{i} - iqr_store.Vbar{i}(end),'c'); hold on; 
    k = k +l+1; 
end
hold on; k = 0; 
for i =  1:length(DDP_store.Vbar) 
    l = length(DDP_store.Vbar{i}); 
    semilogy(k:k+l-1,DDP_store.Vbar{i} - DDP_store.Vbar{i}(end),'r'); hold on; 
    k = k +l+1; 
end






