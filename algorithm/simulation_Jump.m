function simulation_Jump(q,f_s,robot_params, params,save)

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
q0 = q(:,1);
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
height = robot_params.bodyHeight - robot_params.hiplinkLength ...
    -robot_params.kneelinkLength;
% plot([-2 1],[-0.404, -0.404]);
% wall = 0.25; high =  params.hght; 
% plot([-2 wall],[-0.404, -0.404]);
% plot([wall 2],[-0.404 + high, -0.404 + high]);
% plot([wall wall],[-0.404, -0.404 + high]);
plot([-1 2],[-0.404, -0.404])
axis([-1 2 -1 1]);
axis equal
axis manual

%1 is front 2 is back 
y{1} = [f_s(:,1,1);0]; 
y{2} = [f_s(:,1,2);0];

pos{1} =[knee_ft_glob(1,2),knee_ft_glob(2,2),0]';
pos{2} = [knee_bc_glob(1,2),knee_bc_glob(2,2),0]';

for i=1:2
    quad.force(i) = drawArrow(.1, .015,[0 0 0]);
end
1==1;
%% animation
ind = floor(linspace(1,params.len_sim,params.len_sim/10));
q = q(:,ind);
f_s = f_s(:,ind,:); 
for k=1:size(q,2)
    trunk_glob = T_body(q(1:3,k))*trunk_loc;
    hip_ft_glob =  T_body(q(1:3,k))*body_T_hip(robot_params.hip_x(1),q(4,k))*hip_loc;
    knee_ft_glob = T_body(q(1:3,k))*body_T_hip(robot_params.hip_x(1),q(4,k))*hip_T_knee(q(5,k))*knee_loc;
    hip_bc_glob =  T_body(q(1:3,k))*body_T_hip(robot_params.hip_x(2),q(6,k))*hip_loc;
    knee_bc_glob = T_body(q(1:3,k))*body_T_hip(robot_params.hip_x(2),q(6,k))*hip_T_knee(q(7,k))*knee_loc;
    
    h(1).XData = trunk_glob(1,:); h(1).YData = trunk_glob(2,:);
    h(2).XData = hip_ft_glob(1,:); h(2).YData = hip_ft_glob(2,:);
    h(3).XData = knee_ft_glob(1,:); h(3).YData = knee_ft_glob(2,:);
    h(4).XData = hip_bc_glob(1,:); h(4).YData = hip_bc_glob(2,:);
    h(5).XData = knee_bc_glob(1,:); h(5).YData = knee_bc_glob(2,:);
    
    
    pos{1} =[knee_ft_glob(1,2),knee_ft_glob(2,2),0]';
    pos{2} = [knee_bc_glob(1,2),knee_bc_glob(2,2),0]';
    
    y{1} = [f_s(:,k,1);0]; 
    y{2} = [f_s(:,k,2);0];

    for i=1:2 
        updateArrow(quad.force(i),pos{i},y{i}/200);
    end    
    frame{k} = getframe(f);
    im{k} = frame2im(frame{k});
    drawnow;
    pause(0.05)
    
    
end
1==1; 
if save
    for k=1:size(q,2)   
        [imind,cm] = rgb2ind(im{k},256);
        if k == 1
            imwrite(imind,cm,params.filename,'gif', 'Loopcount',inf,'DelayTime',params.dt);
        elseif k==size(q,2)
            imwrite(imind,cm,params.filename,'gif','WriteMode','append','DelayTime',2);
        else
            imwrite(imind,cm,params.filename,'gif','WriteMode','append','DelayTime',params.dt);
        end
    end
end 
end