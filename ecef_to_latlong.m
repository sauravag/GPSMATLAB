%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Convert ECEF co-ordinates to Geodetic        %
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% All i/o units are specified in brackets %
%   Conventions:
%                1. WGS-84 system used for geodetic model
%                2. Geodetic Coordinates are used in the form lat(rad)/long(rad)/alt(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References: 
%    1. J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates to geodetic
%       coordinates," Aerospace and Electronic Systems, IEEE Transactions on, vol. 30, pp. 957–961, 1994.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
%        1. pos_geodetic:  lat(rad)/long(rad)/alt(m)
% Input:
%        2. pos_ecef: ECEF position (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_geodetic] = ecef_to_latlong(pos_ecef);
   
    % WGS84 ellipsoid constants:
    a = 6378137;
    b = 6356752.3142; % semi-minor axis
    e = 8.1819190842622e-2;
    e1 = sqrt(6.73949674228e-3);
    
    x = pos_ecef.x;
    y = pos_ecef.y;
    z = pos_ecef.z;
    
    r = sqrt(x^2+y^2);
    E = sqrt(a^2-b^2);
    F = 54*b^2*z^2;
    G = r^2+(1-e^2)*z^2-e^2*E^2;
    C = (e^4*F*r^2)/G^3;
    S = nthroot(1+C+sqrt(C^2+2*C),3);
    P = F/(3*(S + 1/S + 1)^2*G^2);
    Q = sqrt(1+2*e^4*P);
    r_o = -P*e^2*r/(1+Q)+ sqrt(a^2/2*(1+1/Q)- P*(1-e^2)*z^2/(Q*(1+Q))- P*r^2/2);
    U = sqrt((r-e^2*r_o)^2+z^2);
    V = sqrt((r-e^2*r_o)^2+(1-e^2)*z^2);
    z_o = b^2*z/(a*V);
    altitude = U*(1-b^2/(a*V));
    latitude = atan((z+e1^2*z_o)/r);
    longitude = atan2(y,x);
    
    pos_geodetic = struct('lat',latitude,'long',longitude,'alt',altitude);
     
end