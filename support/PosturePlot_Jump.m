function PosturePlot_Jump(q, robot_params)

%% help functions
RotY = @(t)[cos(t)   sin(t)  0;
            -sin(t)  cos(t)  0
            0        0       1];
            
Trans =@(v)[1        0        v(1);
            0        1        v(2);
            0        0        1];

T_body = @(c)Trans([c(1), c(2)])*RotY(c(3));         % body frame w.r.t. inertia frame

body_T_hip = @(hip_x,t) Trans([hip_x, 0])*RotY(t);  % front knee frame w.r.t. body frame

hip_T_knee = @(t) Trans([0 -robot_params.hiplinkLength])*RotY(t);  
           
%% link data w.r.t. local frame
bodyLength = robot_params.bodyLength+0.04;
bodyHeight = robot_params.bodyHeight+0.04;
hipLength = robot_params.hiplinkLength;
kneeLength = robot_params.kneelinkLength;
linkwidth = 0.03;
% pos of lower left corner of the link
trunk_loc = [-bodyLength/2, -bodyLength/2,  bodyLength/2,  bodyLength/2; 
             bodyHeight/2,  -bodyHeight/2,  -bodyHeight/2,  bodyHeight/2;
             1              1,              1,              1];

hip_loc = [-linkwidth/2,    -linkwidth/2,  linkwidth/2,  linkwidth/2;
           0,               -hipLength,    -hipLength,   0;
           1                1,             1,               1];

knee_loc = [-linkwidth/2,    -linkwidth/2,  linkwidth/2,  linkwidth/2;
           0,               -kneeLength,    -kneeLength,   0;
           1                1,             1,               1];

% plot initial position
q0 = q;
trunk_glob = T_body(q0(1:3))*trunk_loc;
hip_ft_glob =  T_body(q0(1:3))*body_T_hip(robot_params.hip_x(1),q0(4))*hip_loc;
knee_ft_glob = T_body(q0(1:3))*body_T_hip(robot_params.hip_x(1),q0(4))*hip_T_knee(q0(5))*knee_loc;
hip_bc_glob =  T_body(q0(1:3))*body_T_hip(robot_params.hip_x(2),q0(6))*hip_loc;
knee_bc_glob = T_body(q0(1:3))*body_T_hip(robot_params.hip_x(2),q0(6))*hip_T_knee(q0(7))*knee_loc;
f=figure;
h(1)=patch(trunk_glob(1,:),trunk_glob(2,:),'red');
hold on;
h(2)=patch(hip_ft_glob(1,:),hip_ft_glob(2,:),'blue');
h(3)=patch(knee_ft_glob(1,:),knee_ft_glob(2,:),'green');
h(4)=patch(hip_bc_glob(1,:),hip_bc_glob(2,:),'blue');
h(5)=patch(knee_bc_glob(1,:),knee_bc_glob(2,:),'green');
% title(sprintf('current time index %d', index));
height = robot_params.bodyHeight - robot_params.hiplinkLength ...
    -robot_params.kneelinkLength;
plot([-1 2],[-0.404, -0.404]);
% plot([-1 0.45],[-0.404, -0.404]);
% plot([0.45 2],[-0.404*-height, -0.404*-height]);
% plot([0.45 0.45],[-0.404, -0.404 * -height]);
axis([-1 2 -1 1]);
axis equal
axis manual

end