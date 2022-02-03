function [Pos_contact,T_front,T_rear]  = fkin_quadruped_2D(q, robot_params)
% q(1):x  q(2):y  q(3):roty
% q(4):theta1 hip  q(5):theta2 knee (front leg)
% q(6):theta1 hip  q(7):theta2 knee (back leg)       
%% basic homogeneous transformation functions
% RotX = @(t)[1    0       0       0;
%             0    cos(t)  -sin(t) 0;
%             0    sin(t)  cost(t) 0;
%             0    0       0       1];
RotY = @(t)[cos(t)   0   sin(t)  0;
            0        1   0       0;
            -sin(t)  0   cos(t)  0;
            0        0   0       1];
% RotZ = @(t)[cos(t)   -sin(t) 0   0;
%             sin(t)   cos(t)  0   0;
%             0        0       1   0;
%             0        0       0   1];
Trans =@(v)[1        0       0   v(1);
            0        1       0   v(2);
            0        0       1   v(3);
            0        0       0   1];

%% kinematics
T_body = Trans([q(1), 0, q(2)])*RotY(q(3));         % body frame w.r.t. inertia frame

body_T_front = Trans([robot_params.hip_x(1), 0, 0])*RotY(q(4))...   % front knee frame w.r.t. body frame
                *Trans([0 0 -robot_params.hiplinkLength])*RotY(q(5));   
            
body_T_rear = Trans([robot_params.hip_x(2), 0, 0])*RotY(q(6))...    % rear knee frame w.r.t. body frame
                *Trans([0 0 -robot_params.hiplinkLength])*RotY(q(7));
            
knee_Pos_knee = [0, 0, -robot_params.kneelinkLength]';                   % foot pos in knee frame
T_front = T_body*body_T_front;
Pos_front = T_body*body_T_front*[knee_Pos_knee;1];   % front foot pos in inertia frame
T_rear = T_body*body_T_rear;
Pos_rear = T_body*body_T_rear*[knee_Pos_knee;1];     % rear foot pos in inertia frame

Pos_front = Pos_front(1:end-1);    % remove the fourth element
Pos_rear = Pos_rear(1:end-1);     % remove the fourth element

Pos_contact = [Pos_front, Pos_rear];
end
