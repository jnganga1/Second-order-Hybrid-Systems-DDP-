% dt = 0.02
% len = 1/dt;
% T = zbar(:,1);
% tspan_bc = T(1)*(1:len)*dt;
% tspan_f1 = T(2)*(1:len)*dt + T(1)*len*dt;
% tspan_ft = T(3)*(1:len)*dt + (T(1)+T(2))*len*dt;
% tspan_f2 = T(4)*(1:len)*dt + (T(1)+T(2)+T(3))*len*dt;
% tspan = [0,tspan_bc,tspan_f1,tspan_ft,tspan_f2];


dt = 0.001;
len = length(xbar);
tspan = dt*(0:len-1);

figure;
subplot(2,2,1)
plot(tspan,ubar(1,:));
subplot(2,2,2)
plot(tspan,ubar(2,:));
subplot(2,2,3)
plot(tspan,ubar(3,:));
subplot(2,2,4)
plot(tspan,ubar(4,:));
title('Joint torque','interpreter','latex');

figure;
subplot(2,2,1)
plot(tspan,ybar(1,:,1));
xlabel('time s','interpreter','latex');
ylabel('$x$','interpreter','latex');
title('$F_{\rm{fr}}$','interpreter','latex');
subplot(2,2,2)
plot(tspan,ybar(2,:,1));
xlabel('time s','interpreter','latex');
ylabel('$z$','interpreter','latex');
title('$F_{\rm{fr}}$','interpreter','latex');

subplot(2,2,3)
plot(tspan,ybar(1,:,2));
xlabel('time s','interpreter','latex');
ylabel('$x$','interpreter','latex');
title('$F_{\rm{bc}}$','interpreter','latex');

subplot(2,2,4)
plot(tspan,ybar(2,:,2));
xlabel('time s','interpreter','latex');
ylabel('$z$','interpreter','latex');
title('$F_{\rm{bc}}$','interpreter','latex');


% figure
% plot(tspan,ubar(2,:),'-','linewidth',1.5);
% ylabel('torque ($N\cdot m$)','interpreter','latex');
% yyaxis right
% plot(tspan,y1bar(2,:),'-','linewidth',1.5);
% xlabel('time ($s$)','interpreter','latex');
% ylabel('GRF ($N$)','interpreter','latex');




% figure;
% subplot(2,2,1)
% plot(tspan,xbar(11,:));
% subplot(2,2,2)
% plot(tspan,xbar(12,:));
% subplot(2,2,3)
% plot(tspan,xbar(13,:));
% subplot(2,2,4)
% plot(tspan,xbar(14,:));
% title('Joint vel','interpreter','latex');


