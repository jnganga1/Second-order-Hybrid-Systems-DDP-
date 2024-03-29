function robot_params = get_miniCheetah_params()
% parameters are retrived from Cheetah Software (newest version up to 1/11/2020)
    robot_params = struct(...
        'CoM_body',  zeros(3,1),...
        'CoM_hip', [0, 0.016, -0.02]',...
        'CoM_knee',[0, 0, -0.061]',...
        'CoM_abad', [0, 0.036, 0]',...
        ...
        'mass_body', 3.3,...
        'mass_hip' , 0.634,...
        'mass_knee' , 0.064,...
        'mass_abad' , 0.54,...
        ... % link rotational inertia w.r.t. its CoM expressed in its own frame
        'I_body' , [11253, 0, 0; 0, 36203, 0; 0, 0, 42673]*1e-6,...          % trunk rotational inertia
        'I_hip' , [1983, 245, 13; 245, 2103, 1.5; 13, 1.5, 408]*1e-6,...     % hip rotational inertia
        'I_knee' , [245, 0, 0; 0, 248, 0; 0, 0, 6]*1e-6,...                  % knee rotational inertia
        'I_abad' , [381, 58, 0.45; 58, 560, 0.95; 0.45, 0.95, 444]*1e-6,...  % abad rotational inertia
        ...
        'bodyLength' , 0.19 * 2,...
        'bodyWidth' , 0.049 * 2,...
        'bodyHeight' , 0.05 * 2,...
        ...
        'hiplinkLength' , 0.209,...
        'kneelinkLength' , 0.195);
end

