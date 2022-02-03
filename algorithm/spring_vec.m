function v = spring_vec(q, robot_params, ind)
%% Help functions
RotY = @(t)[cos(t)   0   sin(t)  0;
            0        1   0       0;
            -sin(t)  0   cos(t)  0;
            0        0   0       1];
Trans =@(v)[1        0       0   v(1);
            0        1       0   v(2);
            0        0       1   v(3);
            0        0       0   1];

%%
switch ind
    case 1 % bs
        leg = 2;
    case 3 % fs
        leg = 1;
end
pos_foot = fkin_quadruped_2D(q,robot_params);
pos_foot = pos_foot(:,leg);
pos_foot(2) = []; % remove the y element

pos_hip = Trans([q(1), 0, q(2)])*RotY(q(3))*[robot_params.hip_x(leg),0,0,1]';
pos_hip([2,4],:) = []; % remove y element and the fourth element

v = pos_hip - pos_foot;
end

