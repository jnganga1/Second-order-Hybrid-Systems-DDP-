import casadi.*
addpath([pwd '\algorithm']);

NJ = 7; % Degree of freedom
NB = 5; % Number of free bodies
NA = 4; % Number of actuations

% state symbolic variables
q_sym = MX.sym('q',[NJ 1]); %'real'); % x, z, pitch, joint angles
q_dot_sym = MX.sym('qd', [NJ 1]);% 'real');
x_sym = [q_sym; q_dot_sym];
y_sym = MX.sym('yr',[2 1]);%'real'); % contact force (2x1) (2D quadruped)
u_sym = MX.sym('u',[NA 1]); %'real'); % actuations

% set up dynamics and retrieve robot_params
% [params, robot_params,DynFns] = set_system_functions_KKTsolve(x_sym, u_sym, '2D');
[params, robot_params,DynFns] = AllTogether_set_system_functions(x_sym, u_sym, '2D');


%Compare them 
x=rand(NJ*2,1); u = rand(NA,1); bbVal= rand(NJ+2,1); 
muVal=rand(NJ,1); 

all_bc = DynFns.bc_dyn_derivs(x,u);
all_bc_kkt = DynFns.kkt_bc_first(x,u);
Err.err1 = full(all_bc - all_bc_kkt);

all_ft = DynFns.ft_dyn_derivs(x,u);
all_ft_kkt = DynFns.kkt_ft_first(x,u);
Err.err2 = full(all_ft - all_ft_kkt);

all_Imp_bc = DynFns.bc_Im_dyn_derivs(x,u);
all_Imp_kkt_bc = DynFns.kkt_bc_Im_first(x,u);
Err.err3 = full(all_Imp_bc - all_Imp_kkt_bc);


all_Imp_ft = DynFns.ft_Im_dyn_derivs(x,u);
all_Imp_kkt_ft = DynFns.kkt_ft_Im_first(x,u);
Err.err4 = full(all_Imp_ft - all_Imp_kkt_ft);

all_fr_dyn = DynFns.fr_dyn_derivs(x,u);
all_kkt_fr = DynFns.kkt_fr_first(x,u);
Err.err5 = full(all_fr_dyn - all_kkt_fr);


%Front Stance
[a,ft_xx,ft_ux]=DynFns.Exp.Ft_stnce(x,u,bbVal);
[b,ext_xx,ext_ux]=DynFns.ExtMod.Ft_stnce(x,u,bbVal);
[c,tens_xx,tens_ux]=DynFns.Tens.Ft_stnce(x,u);
tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
tens_xx = Tens3byVec(tens_xx,bbVal','pre');
tens_ux = Tens3byVec(tens_ux,bbVal','pre');

% [~,~,John]=DynFns.Tens.John(x,u); John = full(John);
% John = reshape(John,[NJ+2,NJ*2,NA]);
% John = Tens3byVec(John,bbVal','pre');
% John 
% tens_ux

% a - b
% b - c
Err.err6 = full(tens_xx - ext_xx);
Err.err7 = full(tens_xx - ft_xx);
Err.err8 = full(tens_ux - ext_ux'); 
Err.err9 = full(tens_ux - ft_ux');

%Back Stance 
[a,bc_xx,bc_ux]=DynFns.Exp.Bc_stnce(x,u,bbVal);
[b,ext_xx,ext_ux]=DynFns.ExtMod.Bc_stnce(x,u,bbVal);
[c,tens_xx,tens_ux]=DynFns.Tens.Bc_stnce(x,u);
tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
tens_xx = Tens3byVec(tens_xx,bbVal','pre');
tens_ux = Tens3byVec(tens_ux,bbVal','pre');

% a - b
% b - c
Err.err10 = full(tens_xx - ext_xx); 
Err.err11 = full(tens_xx - bc_xx);
Err.err12 = full(tens_ux - ext_ux'); 
Err.err13 = full(tens_ux - bc_ux');

%Ft Impact
[a,bc_xx,bc_ux]=DynFns.Exp.Ft_Imp(x,u,bbVal);
[b,ext_xx,ext_ux]=DynFns.ExtMod.Ft_Imp(x,u,bbVal);
[c,tens_xx,tens_ux]=DynFns.Tens.Ft_Imp(x,u);
tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
tens_xx = Tens3byVec(tens_xx,bbVal','pre');
tens_ux = Tens3byVec(tens_ux,bbVal','pre');

% a - b
% b - c
Err.err14 = full(tens_xx - ext_xx); 
Err.err15 = full(tens_xx - bc_xx);
Err.err16 = full(tens_ux - ext_ux'); 
Err.err17 = full(tens_ux - bc_ux');

%BC Impact
[a,bc_xx,bc_ux]=DynFns.Exp.Bc_Imp(x,u,bbVal);
[b,ext_xx,ext_ux]=DynFns.ExtMod.Bc_Imp(x,u,bbVal);
[c,tens_xx,tens_ux]=DynFns.Tens.Bc_Imp(x,u);
tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
tens_xx = Tens3byVec(tens_xx,bbVal','pre');
tens_ux = Tens3byVec(tens_ux,bbVal','pre');

% a - b
% b - c
Err.err18 = full(tens_xx - ext_xx); 
Err.err19 = full(tens_xx - bc_xx);
Err.err20 = full(tens_ux - ext_ux'); 
Err.err21 = full(tens_ux - bc_ux');

%Fr Dyn
[a,fr_xx,fr_ux]=DynFns.Exp.Fr_Dyn(x,u,muVal);
[b,ext_xx,ext_ux]=DynFns.ExtMod.Fr_Dyn(x,u,muVal);
[c,tens_xx,tens_ux]=DynFns.Tens.Fr_Dyn(x,u);
tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
tens_xx = reshape(tens_xx,[NJ,NJ*2,NJ*2]);
tens_ux = reshape(tens_ux,[NJ,NJ,NA]);
tens_xx = Tens3byVec(tens_xx,muVal','pre');
tens_ux = Tens3byVec(tens_ux,muVal','pre');

% a - b
% b - c
Err.err22 = full(tens_xx - ext_xx); 
Err.err23 = full(tens_xx - fr_xx);
Err.err24 = full(tens_ux - ext_ux'); 
Err.err25 = full(tens_ux - fr_ux');


fn = fieldnames(Err);
for k=1:numel(fn)
    msg = sprintf("Well, Number %i Failed",k);
    assert(norm(Err.(fn{k})) < 1e-6,msg)
end



1==1;
numrepts = 50;
exp_time = []; extMod_time = []; tens_time=[]; first_time =[];
for irepts =1:numrepts
    Exp = 0; ExtMod = 0; Tens = 0; 
    first=0; %basis
    %Compare them 
    x=rand(NJ*2,1); u = rand(NA,1); bbVal= rand(NJ+2,1); 
    muVal=rand(NJ,1); 
    
    start =tic; 
    all_bc = DynFns.bc_dyn_derivs(x,u);
    all_ft = DynFns.ft_dyn_derivs(x,u);
    all_bc_Im = DynFns.bc_Im_dyn_derivs(x,u);
    all_ft_Im = DynFns.ft_Im_dyn_derivs(x,u);
    all_fr = DynFns.fr_dyn_derivs(x,u);
    first = first + toc(start);

    %Front Stance
    start =tic; 
    [~,ft_xx,ft_ux]=DynFns.Exp.Ft_stnce(x,u,bbVal);
    Exp = Exp + toc(start);
    
    start = tic;
    [~,ext_xx,ext_ux]=DynFns.ExtMod.Ft_stnce(x,u,bbVal);
    ExtMod = ExtMod + toc(start);
    
    start = tic;
    [~,tens_xx,tens_ux]=DynFns.Tens.Ft_stnce(x,u);
    tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
    tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
    tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
    tens_xx = Tens3byVec(tens_xx,bbVal','pre');
    tens_ux = Tens3byVec(tens_ux,bbVal','pre');
    Tens = Tens + toc(start);

    
    %Back Stance 
    start = tic;
    [~,bc_xx,bc_ux]=DynFns.Exp.Bc_stnce(x,u,bbVal);
    Exp = Exp + toc(start); 
    
    start = tic;
    [~,ext_xx,ext_ux]=DynFns.ExtMod.Bc_stnce(x,u,bbVal);
    ExtMod = ExtMod + toc(start); 
    
    start = tic; 
    [~,tens_xx,tens_ux]=DynFns.Tens.Bc_stnce(x,u);
    tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
    tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
    tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
    tens_xx = Tens3byVec(tens_xx,bbVal','pre');
    tens_ux = Tens3byVec(tens_ux,bbVal','pre');
    Tens = Tens+toc(start);

    %Ft Impact
    start = tic; 
    [~,bc_xx,bc_ux]=DynFns.Exp.Ft_Imp(x,u,bbVal);
    Exp = Exp + toc(start); 
    
    start = tic; 
    [~,ext_xx,ext_ux]=DynFns.ExtMod.Ft_Imp(x,u,bbVal);
    ExtMod = toc(start); 
    
    start = tic;
    [~,tens_xx,tens_ux]=DynFns.Tens.Ft_Imp(x,u);
    tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
    tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
    tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
    tens_xx = Tens3byVec(tens_xx,bbVal','pre');
    tens_ux = Tens3byVec(tens_ux,bbVal','pre');
    Tens = Tens + toc(start);
    
    %BC Impact
    start = tic; 
    [~,bc_xx,bc_ux]=DynFns.Exp.Bc_Imp(x,u,bbVal);
    Exp = Exp + toc(start);
    
    start = tic;
    [~,ext_xx,ext_ux]=DynFns.ExtMod.Bc_Imp(x,u,bbVal);
    ExtMod = ExtMod + toc(start);
    
    start = tic;
    [~,tens_xx,tens_ux]=DynFns.Tens.Bc_Imp(x,u);
    tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
    tens_xx = reshape(tens_xx,[NJ+2,NJ*2,NJ*2]);
    tens_ux = reshape(tens_ux,[NJ+2,NJ,NA]);
    tens_xx = Tens3byVec(tens_xx,bbVal','pre');
    tens_ux = Tens3byVec(tens_ux,bbVal','pre');
    Tens = Tens + toc(start); 

    %Fr Dyn
    start = tic;
    [~,fr_xx,fr_ux]=DynFns.Exp.Fr_Dyn(x,u,muVal);
    Exp = Exp + toc(start); 
    
    start = tic;
    [~,ext_xx,ext_ux]=DynFns.ExtMod.Fr_Dyn(x,u,muVal);
    ExtMod = ExtMod + toc(start); 
    
    start = tic;
    [~,tens_xx,tens_ux]=DynFns.Tens.Fr_Dyn(x,u);
    tens_xx = full(tens_xx);tens_ux = full(tens_ux); 
    tens_xx = reshape(tens_xx,[NJ,NJ*2,NJ*2]);
    tens_ux = reshape(tens_ux,[NJ,NJ,NA]);
    tens_xx = Tens3byVec(tens_xx,muVal','pre');
    tens_ux = Tens3byVec(tens_ux,muVal','pre');
    Tens = Tens+toc(start);
    
    
    first_time(end+1) = first;
    exp_time(end+1) = Exp; 
    extMod_time(end+1) = ExtMod; 
    tens_time(end+1) = Tens; 
end

1==1;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

DIST = 'normal';
labels = {'First Derivs', 'Explicit', 'Extended RNEA', 'Tensor'};
figure; lw = 4;  nbins =10; light =0.2;
FD=histfit(first_time,nbins,DIST);hold on; 
FD(2).DisplayName = labels{1}; FD(2).LineWidth = lw; 
exp=histfit(exp_time,nbins,DIST); 
exp(2).DisplayName = labels{2}; exp(2).LineWidth = lw; 
extmod = histfit(extMod_time,nbins,DIST); 
extmod(2).DisplayName = labels{3}; extmod(2).LineWidth = lw; 
tens = histfit(tens_time,nbins,DIST);
tens(2).DisplayName = labels{4};  tens(2).LineWidth = lw; 
colored = {'b','m','g','r'};
FD(1).FaceColor = colored{1}; FD(1).FaceAlpha = light; 
FD(2).Color = colored{1};
exp(1).FaceColor = colored{2}; exp(1).FaceAlpha = light;
exp(2).Color = colored{2};
extmod(1).FaceColor = colored{3}; extmod(1).FaceAlpha = light; 
extmod(2).Color = colored{3};
tens(1).FaceColor = colored{4}; tens(1).FaceAlpha = light; 
tens(2).Color = colored{4};
xlabel('Logarithm of Time(s)','Interpreter','latex')
Lgd =legend([FD(2),exp(2),extmod(2),tens(2)]);
Lgd.Location = 'best';Lgd.FontSize  = 25;
set(gca,'FontSize',20);
set(gca, 'XScale', 'log');
grid on; grid minor;
title('Derivatives Evaluation (50 evaluation points)')


figure;
all = [first_time' exp_time' extMod_time' tens_time'];
Tme = boxplot(all,'labels',labels); %'whisker',2, 'widths',0.4,'symbol','');
set(gca, 'YScale', 'log'); 
set(gca,'FontSize',30);
set(Tme,'LineWidth', 2);
L = gca;
L.XAxis.TickLabelInterpreter = 'latex';
grid on; grid minor; 
xtickangle(15);
title('Derivatives Evaluation (50 evaluation points)')
1==1;
