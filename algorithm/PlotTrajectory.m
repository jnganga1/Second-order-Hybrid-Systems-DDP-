function PlotTrajectory(xbar, xf, params)
dt = params.dt;
N = params.len_sim;

figure;
subplot(1,2,1)
plot(0:dt:dt*(N-1),xbar(1,:));
hold on
plot(0:dt:dt*(N-1),xf(1,1)*ones(1,N),'r--');
xlabel('t');
ylabel('x');
subplot(1,2,2)
plot(0:dt:dt*(N-1),xbar(8,:));
hold on
plot(0:dt:dt*(N-1),xf(8,1)*ones(1,N),'r--');
xlabel('t');
ylabel('$\dot{x}$','interpreter','latex');

figure;
subplot(1,2,1)
plot(0:dt:dt*(N-1),xbar(2,:));
hold on
plot(0:dt:dt*(N-1),xf(2,1)*ones(1,N),'r--');
xlabel('t');
ylabel('z');
subplot(1,2,2)
plot(0:dt:dt*(N-1),xbar(9,:));
hold on
plot(0:dt:dt*(N-1),xf(9,1)*ones(1,N),'r--');
xlabel('t');
ylabel('$\dot{z}$','interpreter','latex');

figure;
subplot(1,2,1)
plot(0:dt:dt*(N-1),xbar(3,:));
hold on
plot(0:dt:dt*(N-1),xf(3,1)*ones(1,N),'r--');
xlabel('t');
ylabel('$\theta$','interpreter','latex');
subplot(1,2,2)
plot(0:dt:dt*(N-1),xbar(10,:));
hold on
plot(0:dt:dt*(N-1),xf(10,1)*ones(1,N),'r--');
xlabel('t');
ylabel('$\dot{\theta}$','interpreter','latex');
