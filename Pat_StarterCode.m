%% Run (full) iLQR and DDP (using Trad Regularization)
Itgr = 1; %DON'T CHANGE ME 
OtherTests.Run = false; 
OtherTests.randomInit = false;

%Using iLQR method - no regularization
iLQR = 1; Method = 'none';
iLQR_store  = DDP_2DQuadruped(iLQR,Method,0,Itgr,OtherTests);

%using DDP, regularization: Trad, Method: ExtMod 
iLQR = 0; Method = 2; Reg = 1;
DDP_ExtMod_store  = DDP_2DQuadruped(iLQR,Method,Reg,Itgr,OtherTests);

