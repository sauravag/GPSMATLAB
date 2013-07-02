%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convert from Geodetic to ECEF coordinates     
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vector_ECEF] = ned_to_ecef(pos_body_ecef,vector)

% C_b2n = [cos(theta)*cos(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);-sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];

pos_body_geo = ecef_to_latlong(pos_body_ecef);
lat = pos_body_geo.lat;
long = pos_body_geo.long;
C_e2n = [-sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat);-sin(long) cos(long) 0; -cos(lat)*cos(long) -cos(lat)*sin(long) -sin(lat)];

vector_ECEF =   (C_e2n'*vector')';

end