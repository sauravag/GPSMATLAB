%Calculates the gravitational acceleration given the ECEF position

function gaccel_ned = eval_gaccel_ned(pos_body_geo)

mu = 3.986004415e14; %m^3/s^2
a = 6378137; %m; equatorial radius of earth
% b = 6356752.3142; %m; polar radius of earth
J2 = 0.00108263;

pos_body_ecef = latlong_to_ecef(pos_body_geo);

r1 = pos_body_ecef.x;
r2 = pos_body_ecef.y;
r3 = pos_body_ecef.z;
r = sqrt(r1^2+r2^2+r3^2);
Lg = r3/r; %geocentric latitude

% C_b2n = [cos(theta)*cos(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);cos(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);-sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];

pos_body_geo = ecef_to_latlong(pos_body_ecef);

lat = pos_body_geo.lat;

long = pos_body_geo.long;

C_ECEF2NED = [-sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat);-sin(long) cos(long) 0; -cos(lat)*cos(long) -cos(lat)*sin(long) -sin(lat)];

gaccel_ECEF = -mu/r^2*...
    [(1+1.5*J2*(a/r)^2*(1-5*(Lg)^2))*r1/r;...
        (1+1.5*J2*(a/r)^2*(1-5*(Lg)^2))*r2/r;...
            (1+1.5*J2*(a/r)^2*(3-5*(Lg)^2))*r3/r];
        
gaccel_ned = C_ECEF2NED*gaccel_ECEF;


end
%Ref: Page 41, Aircraft Control and Simulation by Stevens & Lewis