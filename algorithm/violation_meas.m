function height_viol = violation_meas(x,robot_params,ind)
% ind = 2 front impact 
% ind = 4 back impact

contact_ref = fkin_quadruped_2D(zeros(robot_params.q_size,1),robot_params); 
% straight-up position as ref

q = x(1:robot_params.q_size,1);
pos_current = fkin_quadruped_2D(q, robot_params);

k_hat = [0 0 1]';
height = k_hat'*(pos_current - contact_ref); % project distance vecotr to z direction
switch ind         % evaluate the transition ind
    case 2     % front impact. measure front leg
        height_viol = height(1);
    case {4,5}    % back impact. measure back leg
        height_viol = height(2);
end
end

