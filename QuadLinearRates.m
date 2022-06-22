%pdf: 
% http://people.whitman.edu/~hundledr/courses/M467F06/ConvAndError.pdf
%Setting Defaults makes life easier
set(groot,'defaultLineLineWidth',4)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',15)

K = 2;

linear = @(n,k) 1/ (n^k);
quad = @(n) (1.3)^(-2^n + 6);

numRepts = 100;
L = [];
Q = L;  
for i=1:1:numRepts
    L(end+1) = linear(i,K); 
    Q(end+1) = quad(i); 
end

Lchange= L;
LsuLinear = abs( Lchange(2:end) ./ Lchange(1:end-1));
LQuad = abs(Lchange(2:end) ./ Lchange(1:end-1).^2  );


Qchange= Q;
QsuLinear = abs(Qchange(2:end) ./ Qchange(1:end-1));
QQuad = abs(Qchange(2:end) ./ Qchange(1:end-1).^2);
%
names={'Linear','Quad'}; 
colors = {'r','b'};
figure; subplot(2,1,1);
plot(LsuLinear,'DisplayName',names{1},'Color',colors{1});
hold on;
plot(LQuad,'DisplayName',names{2},'Color',colors{2});
legend; grid on; grid minor;
title('Defintion 2, Function: $\frac{1}{n^2}$','Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('Convergence','Interpreter','latex');
ylim([0 2])

subplot(2,1,2);
plot(QsuLinear,'DisplayName',names{1},'Color',colors{1});
hold on;
plot(QQuad,'DisplayName',names{2},'Color',colors{2});
legend; grid on; grid minor;
title('Definition 2, Function: $10^{{-2}^n}$','Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('Convergence','Interpreter','latex');
% ylim([0 30])


%%
figure; 
semilogy(Lchange,'DisplayName',names{1}); hold on
semilogy(Qchange,'DisplayName',names{2});
grid on; grid minor;
legend
xlabel('$n$','Interpreter','latex');
ylabel('Suboptimality','Interpreter','latex');



%%
names={'$\frac{1}{n^2}$','$10^{{-2}^n}$'}; 
figure; subplot(2,1,1); 
semilogy(LsuLinear,'DisplayName',names{1},'Color',colors{1});
hold on;
semilogy(QsuLinear,'DisplayName',names{2},'Color',colors{2});
title('Linear Convergence Plot','Interpreter','latex');
legend; grid on; grid minor;
xlabel('$n$','Interpreter','latex');
ylabel('Suboptimality','Interpreter','latex');

subplot(2,1,2); 
semilogy(LQuad,'DisplayName',names{1},'Color',colors{1});
hold on;
semilogy(QQuad,'DisplayName',names{2},'Color',colors{2});
title('Quadratic Convergence Plot','Interpreter','latex');
legend; grid on; grid minor;
xlabel('$n$','Interpreter','latex');
ylabel('Suboptimality','Interpreter','latex');



