function [contour] = source_sink(frame)

gim = rgb2gray(frame);

main_sink_loc = [round(width(gim)/2), round(length(gim)/3+5)]; 
main_source_loc = [round(width(gim)/2), round(length(gim)*2/3-8)];
mid_ref = [abs(main_sink_loc(1)+main_source_loc(1))/2, abs(main_sink_loc(2)+main_source_loc(2))/2];


main_xdiff = main_source_loc(:,1)-main_sink_loc(:,1);
main_ydiff = main_source_loc(:,2)-main_sink_loc(:,2);
main_angle = atan(main_ydiff/main_xdiff); 
perp_angle = main_angle - 90*pi/180; 
steps = 12; stepsize = 3; 
n = 0; 

% This while loop just iterates and adjusts the central reference points
% the points are named based on their location around the LV (e.g. bottom
% right, bottom_left etc.)
while n <= 3

bottom_rightx = zeros(1,steps); bottom_righty = zeros(1,steps); 
bottom_right = zeros(1,steps);
bottom_rightx(1) = main_source_loc(:,1); 
bottom_righty(1) = main_source_loc(:,2);

bottom_leftx = zeros(1,steps); bottom_lefty = zeros(1,steps); 
bottom_left = zeros(1,steps);
bottom_leftx(1) = main_source_loc(:,1); 
bottom_lefty(1) = main_source_loc(:,2);

bottom_x = zeros(1,steps); bottom_y = zeros(1,steps); 
bottom = zeros(1,steps);
bottom_x(1) = main_source_loc(:,1); 
bottom_y(1) = main_source_loc(:,2);

top_rightx = zeros(1,steps); top_righty = zeros(1,steps); 
top_right = zeros(1,steps);
top_rightx(1) = main_sink_loc(:,1); 
top_righty(1) = main_sink_loc(:,2);

top_leftx = zeros(1,steps); top_lefty = zeros(1,steps); 
top_left = zeros(1,steps);
top_leftx(1) = main_sink_loc(:,1); 
top_lefty(1) = main_sink_loc(:,2);

top_x = zeros(1,steps); top_y = zeros(1,steps); 
top = zeros(1,steps);
top_x(1) = main_sink_loc(:,1); 
top_y(1) = main_sink_loc(:,2);

mid_rightx = zeros(1,steps); mid_righty = zeros(1,steps); 
mid_right = zeros(1,steps);
mid_rightx(1) = mid_ref(:,1); 
mid_righty(1) = mid_ref(:,2);

mid_leftx = zeros(1,steps); mid_lefty = zeros(1,steps); 
mid_left = zeros(1,steps);
mid_leftx(1) = mid_ref(:,1); 
mid_lefty(1) = mid_ref(:,2);

% Set threshold for finding supplemental sink/source locations
thresh = 0.4;

% figure
% imshow(gim)
% hold on
% scatter(main_source_loc(1), main_source_loc(2),'*'); 
% scatter(main_sink_loc(1), main_sink_loc(2),'*');
% scatter(mid_ref(1), mid_ref(2),'x'); 

% Find the average value (3x3 neighborhood) at each point extending 
% from the main reference points within the LV 
for i = 2:steps
    % From source
    bottom_rightx(i) = round(bottom_rightx(i-1) + stepsize*cos(perp_angle)); 
    bottom_righty(i) = round(bottom_righty(i-1) + stepsize*sin(perp_angle));
    x_hood_br = [-2, -1, 0, 1, 2] + bottom_rightx(i); 
    y_hood_br = [-2, -1, 0, 1, 2] + bottom_righty(i);

    bottom_leftx(i) = round(bottom_leftx(i-1) - stepsize*cos(perp_angle)); 
    bottom_lefty(i) = round(bottom_lefty(i-1) - stepsize*sin(perp_angle));
    x_hood_bl = [-2, -1, 0, 1, 2] + bottom_leftx(i); 
    y_hood_bl = [-2, -1, 0, 1, 2] + bottom_lefty(i);

    bottom_x(i) = round(bottom_x(i-1) + stepsize*cos(main_angle)); 
    bottom_y(i) = round(bottom_y(i-1) + stepsize*sin(main_angle));
    x_hood_b = [-2, -1, 0, 1, 2] + bottom_x(i); 
    y_hood_b = [-2, -1, 0, 1, 2] + bottom_y(i);

    % From middle reference point
    mid_leftx(i) = round(mid_leftx(i-1) - stepsize*cos(perp_angle)); 
    mid_lefty(i) = round(mid_lefty(i-1) - stepsize*sin(perp_angle));
    x_hood_ml = [-2, -1, 0, 1, 2] + mid_leftx(i); 
    y_hood_ml = [-2, -1, 0, 1, 2] + mid_lefty(i);

    mid_rightx(i) = round(mid_rightx(i-1) + stepsize*cos(perp_angle)); 
    mid_righty(i) = round(mid_righty(i-1) + stepsize*sin(perp_angle));
    x_hood_mr = [-2, -1, 0, 1, 2] + mid_rightx(i); 
    y_hood_mr = [-2, -1, 0, 1, 2] + mid_righty(i);

    % From sink
    top_rightx(i) = round(top_rightx(i-1) + stepsize*cos(perp_angle)); 
    top_righty(i) = round(top_righty(i-1) + stepsize*sin(perp_angle));
    x_hood_tr = [-2, -1, 0, 1, 2] + top_rightx(i); 
    y_hood_tr = [-2, -1, 0, 1, 2] + top_righty(i);

    top_leftx(i) = round(top_leftx(i-1) - stepsize*cos(perp_angle)); 
    top_lefty(i) = round(top_lefty(i-1) - stepsize*sin(perp_angle));
    x_hood_tl = [-2, -1, 0, 1, 2] + top_leftx(i); 
    y_hood_tl = [-2, -1, 0, 1, 2] + top_lefty(i);

    top_x(i) = round(top_x(i-1) - stepsize*cos(main_angle)); 
    top_y(i) = round(top_y(i-1) - stepsize*sin(main_angle));
    x_hood_t = [-2, -1, 0, 1, 2] + top_x(i); 
    y_hood_t = [-2, -1, 0, 1, 2] + top_y(i);


    for j = 1:5
        for k = 1:5
            hood_vals_br(j,k) = gim(y_hood_br(j), x_hood_br(k)); 
            hood_vals_bl(j,k) = gim(y_hood_bl(j), x_hood_bl(k));
            hood_vals_b(j,k) = gim(y_hood_b(j), x_hood_b(k));

            hood_vals_mr(j,k) = gim(y_hood_mr(j), x_hood_mr(k)); 
            hood_vals_ml(j,k) = gim(y_hood_ml(j), x_hood_ml(k));

            hood_vals_tr(j,k) = gim(y_hood_tr(j), x_hood_tr(k)); 
            hood_vals_tl(j,k) = gim(y_hood_tl(j), x_hood_tl(k));
            hood_vals_t(j,k) = gim(y_hood_t(j), x_hood_t(k));
        end
    end
    bottom_right(i) = mean(mean(hood_vals_br)); 
    bottom_left(i) = mean(mean(hood_vals_bl));
    bottom(i) = mean(mean(hood_vals_b));

    middle_right(i) = mean(mean(hood_vals_mr)); 
    middle_left(i) = mean(mean(hood_vals_ml));

    top_right(i) = mean(mean(hood_vals_tr)); 
    top_left(i) = mean(mean(hood_vals_tl));
    top(i) = mean(mean(hood_vals_t));

end

% Normalizing
bottom_right = bottom_right/max(bottom_right); 
bottom_left = bottom_left/max(bottom_left);
bottom = bottom/max(bottom); 

middle_right = middle_right/max(middle_right); 
middle_left = middle_left/max(middle_left);

top_right = top_right/max(top_right); 
top_left = top_left/max(top_left); 
top = top/max(top); 

% Choose the next point out after passing the threshold to set the
% supplmental sources (i+1) & the current point for sinks (i)
for i = 2:steps
    if bottom_right(i) >= thresh && bottom_right(i-1) <= thresh
        if i == steps
            source_br = [bottom_rightx(i), bottom_righty(i)];
        else
            source_br = [bottom_rightx(i), bottom_righty(i)];
        end
        bottom_right = zeros(1,steps);
        s1 = i; 
    end
    if bottom_left(i) >= thresh && bottom_left(i-1) <= thresh
        if i == steps
            source_bl = [bottom_leftx(i), bottom_lefty(i)];
        else
            source_bl = [bottom_leftx(i), bottom_lefty(i)]; 
        end
        bottom_left = zeros(1,steps);
        s2 = i; 
    end
    if bottom(i) >= thresh && bottom(i-1) <= thresh
        source_b = [bottom_x(i), bottom_y(i)]; 
        bottom = zeros(1,steps);
    end
    
    if middle_right(i) >= thresh && middle_right(i-1) <= thresh
        sink_mr = [mid_rightx(i), mid_righty(i)]; 
        middle_right = zeros(1,steps);
        s3 = i; 
    end
    if middle_left(i) >= thresh && middle_left(i-1) <= thresh
        sink_ml = [mid_leftx(i), mid_lefty(i)]; 
        middle_left = zeros(1,steps);
        s4 = i; 
    end
    

    if top_right(i) >= thresh && top_right(i-1) <= thresh
        sink_tr = [top_rightx(i), top_righty(i)]; 
        top_right = zeros(1,steps);
        s5 = i; 
    end
    if top_left(i) >= thresh && top_left(i-1) <= thresh
        sink_tl = [top_leftx(i), top_lefty(i)]; 
        top_left = zeros(1,steps);
        s6 = i; 
    end
    if max(top) <= 100 && max(top) >= 0
        sink_t = [main_sink_loc(1)-abs(main_sink_loc(1)-main_source_loc(1))*2/3, ...
            main_sink_loc(2)-abs(main_sink_loc(2)-main_source_loc(2))*2/3]; 
        top = zeros(1,steps);
    end
    if top(i) >= thresh && top(i-1) <= thresh
        sink_t = [top_x(i), top_y(i)]; 
        top = zeros(1,steps);
    end
end

if n <= 3
    % Reposition main sink/source between identified points
    main_sink_loc(1) = sink_tl(1) + abs(sink_tl(1)-sink_tr(1))*4/8;
    main_sink_loc(2) = sink_tl(2) - abs(sink_tl(2)-sink_tr(2))*4/8;
    main_source_loc(1) = source_bl(1) + abs(source_bl(1)-source_br(1))*4/8;
    main_source_loc(2) = source_bl(2) - abs(source_bl(2)-source_br(2))*4/8;
    main_xdiff = main_source_loc(:,1)-main_sink_loc(:,1);
    main_ydiff = main_source_loc(:,2)-main_sink_loc(:,2);
    main_angle = atan(main_ydiff/main_xdiff); 
    perp_angle = main_angle - 90*pi/180; 

    bottom_dist = sqrt((source_b(1)-main_source_loc(1))^2+(source_b(2)-main_source_loc(2))^2);
    if bottom_dist <= 5

        main_source_loc(1) = main_source_loc(1) - 5*cos(main_angle); 
        main_source_loc(2) = main_source_loc(2) - 5*sin(main_angle); 

    end

    mid_ref = [abs(main_sink_loc(1)+main_source_loc(1))/2, abs(main_sink_loc(2)+main_source_loc(2))/2];
end

% Summing to control number of loops (can be adjusted to see how it changes
% accuracy
n = n+1; 

end

% figure
% imshow(gim)
% hold on
% scatter(main_source_loc(1), main_source_loc(2),'blue', '*'); 
% scatter(main_sink_loc(1), main_sink_loc(2),'red','*');
% scatter(source_bl(1), source_bl(2), 'green', 'x')
% scatter(source_br(1), source_br(2), 'green', 'x')
% scatter(source_b(1), source_b(2), 'green', 'x')
% scatter(sink_tl(1), sink_tl(2), 'green', 'x')
% scatter(sink_tr(1), sink_tr(2), 'green', 'x')
% scatter(sink_t(1), sink_t(2), 'green', 'x')
% scatter(mid_ref(1), mid_ref(2),'red','*');
% scatter(sink_ml(1), sink_ml(2), 'green', 'x')
% scatter(sink_mr(1), sink_mr(2), 'green', 'x')

%end

%% Finding Streamline using superposition velocity

% Setting the strengths of the sinks and sources Q-source, K-sink
Q = 80; %K = 120; Q_sup = -60; K_sup1 = 40; K_sup2 = 100;
eud_thresh = 2; 

for s = 1:2

if s == 1
    p = 1; 
else
    p = -1; 
end

%current_point(1,:) = [length(gim)-(38+s), width(gim)];
current_point(1,:) = [source_b(1)+p*10*sin(main_angle), source_b(2)-p*10*cos(main_angle)];
V_inf = 0.9*Q/sqrt((source_b(1)-main_source_loc(1))^2+(source_b(2)-main_source_loc(2))^2); 
V_infx = V_inf*cos(main_angle+pi); V_infy = V_inf*sin(main_angle+pi);
Sink_xdiff = current_point(1) - main_sink_loc(1); signx5 =  Sink_xdiff/abs(Sink_xdiff);
Sink_ydiff = current_point(2) - main_sink_loc(2); signy5 = Sink_ydiff/abs(Sink_ydiff);
Sink_ang = abs(atan(Sink_ydiff/Sink_xdiff));
eud_Sink = sqrt(Sink_xdiff^2 + Sink_ydiff^2);
n = 1; 

while eud_Sink >= 4 && n <= 40
Q = 80; K = 120; Q_sup = -10; K_sup1 = 10; K_sup2 = 15;

% Supplemental Source 1 (bottom right)
Source1_xdiff = current_point(n,1) - source_br(1); signx1 =  Source1_xdiff/(abs(Source1_xdiff)+1E-10);
Source1_ydiff = current_point(n,2) - source_br(2); signy1 = Source1_ydiff/(abs(Source1_ydiff)+1E-10);
Source1_ang = abs(atan(Source1_ydiff/Source1_xdiff));%-90*pi/180;
eud_Source1 = sqrt(Source1_xdiff^2 + Source1_ydiff^2); 

if eud_Source1 <= eud_thresh
    Q_sup = 0; 
end

Source1_v = Q_sup/(eud_Source1); 
Source1_vx = signx1*Source1_v*cos(Source1_ang); Source1_vy = signy1*Source1_v*sin(Source1_ang);

% Supplmenetal Sink 1 (top right)
Sink1_xdiff = current_point(n,1) - sink_tr(1); signx2 =  -Sink1_xdiff/(abs(Sink1_xdiff)+1E-10);
Sink1_ydiff = current_point(n,2) - sink_tr(2); signy2 = -Sink1_ydiff/(abs(Sink1_ydiff)+1E-10);
Sink1_ang = abs(atan(Sink1_ydiff/Sink1_xdiff));%-90*pi/180;
eud_Sink1 = sqrt(Sink1_xdiff^2 + Sink1_ydiff^2); 

if eud_Sink1 <= eud_thresh
    K_sup1 = 0; 
end

Sink1_v = K_sup1/(eud_Sink1); 
Sink1_vx = signx2*Sink1_v*cos(Sink1_ang); Sink1_vy = signy2*Sink1_v*sin(Sink1_ang);

% Supplmenetal Sink 2 (top left)
Sink2_xdiff = current_point(n,1) - sink_tl(1); signx3 =  -Sink2_xdiff/(abs(Sink2_xdiff)+1E-10);
Sink2_ydiff = current_point(n,2) - sink_tl(2); signy3 = -Sink2_ydiff/(abs(Sink2_ydiff)+1E-10);
Sink2_ang = abs(atan(Sink2_ydiff/Sink2_xdiff));%+90*pi/180;
eud_Sink2 = sqrt(Sink2_xdiff^2 + Sink2_ydiff^2); 

if eud_Sink2 <= eud_thresh
    K_sup1 = 0; 
end

Sink2_v = 3/2*K_sup1/(eud_Sink2); 
Sink2_vx = signx3*Sink2_v*cos(Sink2_ang); Sink2_vy = signy3*Sink2_v*sin(Sink2_ang);

% Supplemental Sink 3 (mid left)
Sink3_xdiff = current_point(n,1) - sink_ml(1); signxm1 =  -Sink3_xdiff/(abs(Sink3_xdiff)+1E-10);
Sink3_ydiff = current_point(n,2) - sink_ml(2); signym1 = -Sink3_ydiff/(abs(Sink3_ydiff)+1E-10);
Sink3_ang = abs(atan(Sink3_ydiff/Sink3_xdiff));%+90*pi/180;
eud_Sink3 = sqrt(Sink3_xdiff^2 + Sink3_ydiff^2); 

if eud_Sink3 <= eud_thresh
    K_sup2 = 0; 
end

Sink3_v = 3/2*K_sup2/(eud_Sink3); 
Sink3_vx = signxm1*Sink3_v*cos(Sink3_ang); Sink3_vy = signym1*Sink3_v*sin(Sink3_ang);

% Supplemental Sink 4 (mid right)
Sink4_xdiff = current_point(n,1) - sink_mr(1); signxm2 =  -Sink4_xdiff/(abs(Sink4_xdiff)+1E-10);
Sink4_ydiff = current_point(n,2) - sink_mr(2); signym2 = -Sink4_ydiff/(abs(Sink4_ydiff)+1E-10);
Sink4_ang = abs(atan(Sink4_ydiff/Sink4_xdiff));%+90*pi/180;
eud_Sink4 = sqrt(Sink4_xdiff^2 + Sink4_ydiff^2); 

if eud_Sink4 <= eud_thresh
    K_sup2 = 0; 
end

Sink4_v = K_sup2/(eud_Sink4); 
Sink4_vx = signxm2*Sink4_v*cos(Sink4_ang); Sink4_vy = signym2*Sink4_v*sin(Sink4_ang);

% Supplemental Source 2 (bottom left)
Source2_xdiff = current_point(n,1) - source_bl(1); signx4 =  Source2_xdiff/(abs(Source2_xdiff)+1E-10);
Source2_ydiff = current_point(n,2) - source_bl(2); signy4 = Source2_ydiff/(abs(Source2_ydiff)+1E-10);
Source2_ang = abs(atan(Source2_ydiff/Source2_xdiff));%+90*pi/180;
eud_Source2 = sqrt(Source2_xdiff^2 + Source2_ydiff^2); 

if eud_Source2 <= eud_thresh
    Q_sup = 0; 
end

Source2_v = 3/2*Q_sup/(eud_Source2); 
Source2_vx = signx4*Source2_v*cos(Source2_ang); Source2_vy = signy4*Source2_v*sin(Source2_ang);

% Main Source 
Source_xdiff = current_point(n,1) - main_source_loc(1); signx5 =  Source_xdiff/(abs(Source_xdiff)+1E-10);
Source_ydiff = current_point(n,2) - main_source_loc(2); signy5 = Source_ydiff/(abs(Source_ydiff)+1E-10);
Source_ang = abs(atan(Source_ydiff/Source_xdiff));
eud_Source = sqrt(Source_xdiff^2 + Source_ydiff^2);
Source_v = Q/(eud_Source); 
Source_vx = signx5*Source_v*cos(Source_ang); Source_vy = signy5*Source_v*sin(Source_ang);

% Main Sink
Sink_xdiff = current_point(n,1) - sink_t(1); signx6 =  -Sink_xdiff/abs(Sink_xdiff);
Sink_ydiff = current_point(n,2) - sink_t(2); signy6 = -Sink_ydiff/abs(Sink_ydiff);
Sink_ang = atan(Sink_ydiff/Sink_xdiff);
eud_Sink = sqrt(Sink_xdiff^2 + Sink_ydiff^2); 
Sink_v = K/(eud_Sink); 
Sink_vx = signx6*Sink_v*cos(Sink_ang); Sink_vy = signx6*Sink_v*sin(Sink_ang);

n = n+1;
% Summing all x and y velocity components
X_vel_comp = [Source1_vx, Sink1_vx, Sink2_vx, Sink3_vx, Sink4_vx, Source2_vx, Source_vx, Sink_vx, V_infx]; 
X_vel = sum(X_vel_comp(:)); X_velsign = X_vel/abs(X_vel);
Y_vel_comp = [Source1_vy, Sink1_vy, Sink2_vy, Sink3_vy, Sink4_vy, Source2_vy, Source_vy, Sink_vy, V_infy];
Y_vel = sum(Y_vel_comp(:)); Y_velsign = Y_vel/abs(Y_vel);
Vel_angle = abs(atan(Y_vel/X_vel)); 
current_point(n,:) = [current_point(n-1,1)+2*X_velsign*cos(Vel_angle), current_point(n-1,2)+2*Y_velsign*sin(Vel_angle)]; 

end

%scatter(current_point(:,1), current_point(:,2), 'blue','*')

if s == 1
    contour = current_point(:,:); 
else
    contour = [contour; flip(current_point(:,:))]; 
end

% scatter(contour(:,1), contour(:,2), 'blue','*')

end