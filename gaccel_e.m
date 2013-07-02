%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the gravitational acceleration given the ECEF position        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
%  Page 41, Aircraft Control and Simulation by Stevens & Lewis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gaccel_body = eval_gaccel_body(pos_body_ecef, pos_body_geo, phi,theta,psi)

mu = 3.986004415e14; %m^3/s^2
a = 6378137; %m; equatorial radius of earth
% b = 6356752.3142; %m; polar radius of earth
J2 = 0.00108263;

r1 = pos_body_ecef(1);
r2 = pos_body_ecef(2);
r3 = pos_body_ecef(3);
r = sqrt(r1^2+r2^2+r3^2);
Lg = r3/r; %geocentric latitude

gaccel_ECEF = -mu/r^2*...
    [(1+1.5*J2*(a/r)^2*(1-5*(Lg)^2))*r1/r;...
        (1+1.5*J2*(a/r)^2*(1-5*(Lg)^2))*r2/r;...
            (1+1.5*J2*(a/r)^2*(3-5*(Lg)^2))*r3/r];
end
