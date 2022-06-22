function [params, callback_params] = set_sim_params(params)
%% simulation parameters
dt = .001;  
seq = ["bs","f2","fs","f2","bs","f1","fs"];%,"f2"]; % Specify mode sequence
seq = ["f2","bs","f1","fs","f2","bs"];%"f1"];%,"fs"];%"f2","bs"]; % Specify mode sequence
num_mode = length(seq);     % num of modes in total
num_mode_gait = 4;          % num of modes in one gait

P = zeros(1,num_mode);   % encode mode seq with integer
T_mode = zeros(1,num_mode);   % assign time interval for each mode
len_mode = zeros(1,num_mode);
for i = 1:num_mode
    switch seq(i)
        case "bs"
            P(i) = 1;
            T_mode(i) = 0.072;
        case "f1"
            P(i) = 2;
            T_mode(i) = 0.050;
        case "fs"
            P(i) = 3;
            T_mode(i) = 0.0720;
        case "f2"
            P(i) = 4;
            T_mode(i) = 0.050;
    end
    len_mode(i) = T_mode(i)/dt;
end
T_mode(1) = dt;
len_mode(1) = 1;           % for initial condition
len_sim = sum(len_mode);   % overall simulation length 

params.seq = seq;
params.dt = dt;
params.P = P;       % mode index set
params.len_sim = len_sim;
params.num_mode_gait = num_mode_gait;
params.num_mode= num_mode;
params.num_tdconstr = get_num_tdconstr(seq);
params.len_mode = len_mode;

params.beta = .5;
params.gamma = .01;
params.Debug = 1;
params.OK = true; % DDP success indicator (regularization and line search)

callback_params.plot_V_prediction = 1;
callback_params.pause = 1;
callback_params.UserAdvanceIteration = 0;
end

function num_tdconstr = get_num_tdconstr(seq)
num_tdconstr = nnz(strcmp(seq,"f2"))+nnz(strcmp(seq,"f1"));
end
