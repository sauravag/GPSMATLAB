%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pioneer UAV Flight Mechanics Model        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_ac_ecef, pos_ac_geodetic,Acceleration_body,AngV,X] = FlightMech_Model(dt,X,dth,de,da,dr)

if (de > 20)
de = 20;
elseif (de<-20)
de = -20;
end
if (da > 20)
da = 20;
elseif (da<-20)
da = -20;
end
if (dr > 20)
dr = 20;
elseif (dr<-20)
dr = -20;
end

% Calculating new state by integrating
[X,DX,Acceleration_body,AngV] = integ_RK4(dt,X,dth,de,da,dr);

% Geodetic Position
pos_ac_geodetic = struct('lat',X(10),'long',X(11),'alt',X(12)); 

% ECEF position
pos_ac_ecef =latlong_to_ecef(pos_ac_geodetic); 

% Velocity in ECEF Coordinates
Velocity_ecef = ned_to_ecef(pos_ac_ecef,[X(1) X(2) X(3)]);

end



