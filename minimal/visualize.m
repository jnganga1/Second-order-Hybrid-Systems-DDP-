function  visualize( hybrid_trajectory, model)
%visualize visualize rigid body simulation model


    bigAnim = figure(198)
    clf
    %f = subplot(211)
    hold on
    % set(f,'Name','3D-Output of a bounding quadruped');  % Window title
    set(bigAnim,'Color','w');         % Background colo
    set(bigAnim,'Renderer','OpenGL')
    drawFloor()

    l = model.p_hip{1}(1) - model.p_hip{4}(1);
    w = -model.p_hip{1}(2) + model.p_hip{4}(2);
    h = w/3;

    quad.body = drawBox([l w h],-model.p_hip{3}',[225,124,22]/256*1.2);
    for i = 1:4
        quad.hip(i) = drawCylinder(.015, model.l1,[92 51 23]/256*1.5);
        quad.shank(i) = drawCylinder(.015, model.l2,[92 51 23]/256*1.5);
        quad.force(i) = drawArrow(.1 , .02 ,[1 0 0]);
    end

    camproj('perspective');
    camtarget([0, 0, +0.6])
    view([-85 8])
    %campos ([5, -10, +2.6]);
    camup ([0,0,1]);
    camva(.7)
    axis off
    box off
    axis equal
    light('Position',[5 20 12 ],'Style','local','Color',[.7 .7 .7]);
    drawnow();

    updateObject(quad.body,[0 0 1]',eye(3));


    Rbase = expm(cross([0 pi/2 0]));


    for k = 1:length(hybrid_trajectory)
       control_params = hybrid_trajectory{k}.control_params;
       x_c           = hybrid_trajectory{k}.x_c; % continuously updating part of state
       x_d           = hybrid_trajectory{k}.x_d; % discretely updating part of state
       t             = hybrid_trajectory{k}.t;


       size_x_c = size(x_c);

       for n = 1:size_x_c(1)
          p    = x_c(n,5:7)';
          quat = x_c(n,1:4)';
          R    = quatToR(quat);
          [f_s]   = control_params.controller(t(n), x_c(n,:)', x_d, model, control_params);

          updateObject(quad.body,p, R);

          for i = 1:4
              %body_pts(2*i  ,:) = pHip'; 
              if x_d.leg_state{i} == 0
                hideObject(quad.hip(i));
                hideObject(quad.shank(i));
                hideObject(quad.force(i));

              else
               
                [~, qn] = getTorque(x_c(n,:)', x_d, model, i, [0 0 0]');
                [pFoot, Rshank, pKnee,RHip] = fwdKin(qn,model);
               
                updateObject(quad.hip(i),p+R*model.p_hip{i},R*RHip*Rbase);
                updateObject(quad.shank(i),p+R*model.p_hip{i}+R*pKnee,R*Rshank*Rbase);

                updateArrow(quad.force(i),p+R*(model.p_hip{i}+pFoot),f_s{i}/500);

              end
          end
          figure(198)
          camtarget([p(1:2)' .6]);
          drawnow;
       end
    end

end