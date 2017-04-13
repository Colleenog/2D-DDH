% function bisector = bisector(x,y,radius,angle)
% 
% % x=center_angle_x ; 
% % y=center_angle_y ;
% % radius = radius_angle; 
% % angle = alpha_tmp; angle = angle_b;
% angle_new = angle - 90;
% angle = angle * pi /180;
% m1 = -tan(angle-pi/2);
% x1 = x + radius/2 * sin(angle);
% y1 = y + radius/2 * cos(angle);
% %angle_new = angle_new .* pi/180;
% %m = -tan(angle_new-pi/2);
% m = -1/m1;
% c =  (y1-m*x1);
% bisector = [m,c];
% end

function [bisector1, bisector2, center_FH, radius_FH] = ...
    bisector (x,y,radius,angle_a,angle_b)


% x=center_angle_x ; 
% y=center_angle_y ;
% radius = radius_angle; 
% angle_a = alpha_tmp; angle_b = angle_b;

angle = angle_a * pi /180;
m1 = -tan(angle-pi/2);
x1 = x + radius/2 * sin(angle);
y1 = y + radius/2 * cos(angle);
m5 = -1/m1;
c5 =  (y1-m5*x1);
bisector1 = [m5,c5];

angle = angle_b * pi /180;
m1 = -tan(angle-pi/2);
x1 = x + radius/2 * sin(angle);
y1 = y + radius/2 * cos(angle);
m6 = -1/m1;
c6 =  (y1-m6*x1);
bisector2 = [m6,c6];

center_x = (c6-c5)/(m5-m6);
center_y = m5 * center_x + c5;
center_FH = [center_x,center_y];

radius_FH = ( (x1-center_x)^2+(y1-center_y)^2 )^(1/2);
end