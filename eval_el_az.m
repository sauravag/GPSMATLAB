%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate Satellite Azimuth and Elevation from User to Satellite       
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
% Montenbruck and Gill pg. 37 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elevation_sv, azimuth_sv, is_visible] = eval_el_az(user_pos_geodetic,user_pos_ecef,sat_pos_ecef);

user_lat = user_pos_geodetic.lat;
user_long = user_pos_geodetic.long;

e_E = [-sin(user_long) cos(user_long) 0];
e_N = [-sin(user_lat)*cos(user_long) -sin(user_lat)*sin(user_long) cos(user_lat)];
e_Z = [cos(user_lat)*cos(user_long) cos(user_lat)*sin(user_long) sin(user_lat)];

E = [e_E;e_N;e_Z];

s = E * [(sat_pos_ecef.x-user_pos_ecef.x); (sat_pos_ecef.y - user_pos_ecef.y); (sat_pos_ecef.z - user_pos_ecef.z)];

azimuth_sv = atan2(s(1),s(2));
elevation_sv = atan2(s(3),sqrt(s(1)^2+s(2)^2));

% check if satellite is visible from the user

is_visible= 0; % default assumed satellite not visible

if (0.17453< elevation_sv && elevation_sv < (3.14159 - 0.17453))
         is_visible = 1; % incase it is visible, value gets changed to 1
end;


end