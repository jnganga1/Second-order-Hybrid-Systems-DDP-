%% Add noise to ubar - Run iLQR/DDP for 20 Iterations 
close all; clear; clc; 
% Add noise to ubar -- some weiner noise.
% Run iLQR/DDP for 50 Iterations and get time - makes figure A 
Itgr = 1;
OtherTests.Run = true; 
OtherTests.randomInit = false;
OtherTests.maxItr = 100;
OtherTests.dt = 1e-3;

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
    OtherTests.scalarNse = 5;


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

%     iLQR_Box{idx,i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %iLQR

    iLQR = 0; Method = 2; Reg = 1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 

%     ExtMod_Box{idx,i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP,ExtMod

    iLQR = 0; Method = 1; Reg = 1;
    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 

%     Exp_Box{idx,i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Exp
    iLQR = 0; Method = 3; Reg = 1;

    OtherTests.iLQR = iLQR; OtherTests.Method = Method; 
    OtherTests.RegMethod = Reg; 
    Tens_Box{i} = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests); %DDP, Tens
end 

% figure; 
% plot(iLQR_Box{1}.Vbar{:},'DisplayName','iLQR'); hold on; 
% plot(ExtMod_Box{1}.Vbar{:},'DisplayName','ExtMod');
% plot(Exp_Box{1}.Vbar{:},'DisplayName','Exp');
% plot(Tens_Box{1}.Vbar{:},'DisplayName','Tens');
% legend

save('MatData/VaryU_Noise_ShortRuns_Tens_fixed_ReModded2.mat');
