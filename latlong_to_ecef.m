%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convert from Geodetic to ECEF coordinates     
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos_ecef] = latlong_to_ecef(user_pos_geodetic);

a = 6378137.0;% semi major axis in m

e = sqrt(6.69437999014e-3);

lat = user_pos_geodetic.lat;

long = user_pos_geodetic.long;

alt = user_pos_geodetic.alt;

X = sqrt(1 - e^2*(sin(lat))^2);

N = a/X;

x = (N + alt)*cos(lat)*cos(long);

y = (N + alt)*cos(lat)*sin(long);

z = (N*(1-e^2)+ alt)*sin(lat);

pos_ecef = struct('x',x,'y',y,'z',z);

end
