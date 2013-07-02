%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pioneer UAV Flight Mechanics Model        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
%            1. "Aircraft Control and Simulation", B.L Stevenson & Frank M. Lewis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:     
%        1. DX: array containing derivative of state vector
%        2. Acceleration_body: acceleration along body axis
%        3. AngV : Angular velocity of aircraft about body axis
% Inputs:
%        1. X: state vector
%        2. dth: throttle (95-500)
%        3. de: elevator deflection (degrees)
%        4. da: aileron deflection (degrees)
%        5. dr: rudder deflection (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[DX,Acceleration_body,AngV] = eval_state_derivative(X,dth,de,da,dr)

%** ADDITIONAL PARAMTERS FOR CALCULATION **%
rtd = 57.29577951; % radian to degree conversion factor
dtr = 0.017453293; % degree to radian conversion factor
Pi =  3.141592654;
Re = 6378137;          %   WGS-84 equatorial radius (m).  
OMEGA_dot_e = 7.292115e-5;%earth's rotation rate (rad/sec)
f = 1/298.257223563;        %   WGS-84 Flattening.
e = sqrt(f*(2 - f));        %   Eccentricity.

pos_body_geo = struct('lat',X(10),'long',X(11),'alt',X(12));

C_b2n = [cos(X(5))*cos(X(6)) -cos(X(4))*sin(X(6))+sin(X(4))*sin(X(5))*cos(X(6)) sin(X(4))*sin(X(6))+cos(X(4))*sin(X(5))*cos(X(6));cos(X(5))*sin(X(6)) cos(X(4))*cos(X(6))+sin(X(4))*sin(X(5))*sin(X(6)) -sin(X(4))*cos(X(6))+cos(X(4))*sin(X(5))*sin(X(6));-sin(X(5)) sin(X(4))*cos(X(5)) cos(X(4))*cos(X(5))];

R_N = Re*(1-e^2)/(1-e^2*(sin(X(10)))^2)^1.5;

R_E = Re/(1-e^2*(sin(X(10)))^2)^0.5; 

w_ien = OMEGA_dot_e*[cos(X(10));0;-sin(X(10))];

w_enn = [X(2)/(R_E+X(12));-X(1)/(R_N+X(12));-X(2)*tan(X(10))/(R_E+X(12))];

% gravitational acceleration
gaccel_ned = eval_gaccel_ned(pos_body_geo);
g_ln = gaccel_ned - OMEGA_dot_e^2*((Re+X(12))/2)*[sin(2*X(10));0;1+cos(2*X(10))];


% Aircraft Parameters 
mass = 190;% [kg]
Ixx = 47.22;% [kg-m^2]
Iyy = 90.84;% [kg-m^2]         
Izz = 111.48;% [kg-m^2]   
Ixz = -6.64;% [kg-m^2]         
Sref = 2.826; % Reference Surface Area [m^2]   
bref = 5.15;% Reference Span  [m]       
cref = 0.55;% Reference Chord Length [m]        
xcg = 1; %?? CG location relative to wing leading edge, expressed as a fraction of aerodynamic chord length        
heng = 0; % ??? Engine angular momentum, assumed fixed 


% Easy representation of moments of interia
Gamma = Ixx * Izz - (Ixz  * Ixz);
C1 = ((Iyy - Izz)  * Izz  - (Ixz * Ixz))/ Gamma;
C2 = ((Ixx - Iyy + Izz ) * Ixz ) / Gamma;
C3 = Izz / Gamma;
C4 = Ixz / Gamma;
C5 = (Izz - Ixx) / Iyy;
C6 = Ixz / Iyy ;
C7 = 1 / Iyy;
C8 = (Ixx * (Ixx - Iyy ) + Ixz * Ixz) / Gamma;
C9 = Ixx / Gamma;

Vbody = C_b2n'*[X(1);X(2);X(3)]; % convert velocity from ned to body frame to calculate aerodynamic forces

Vt = norm(Vbody);
BETA = rtd * asin(Vbody(2)/Vt);
ALPHA = rtd * atan(Vbody(3)/Vbody(1));
DX = linspace(0,0,12);
[qbar,rho] = atmosphere(X(12),Vt);
Thrust = eval_Thrust(rho,Vt,dth);
[CD_fa,CD_fade,CL_fa,CL_fade,CY_fbdr,Cl_fada,Cm_fa,Cm_fade,Cn_fbdr,Cn_fada] = eval_aerod_coeff(de,da,dr,ALPHA,BETA);

%Total force coefficients 

% Cx_tot
CX_tot = CD_fa + CD_fade;

% Cy_tot
CY_tot = CY_fbdr;

%/ Cz_tot
CZ_tot = CL_fa + CL_fade;


% Total moment coefficients

% Cl_tot
Cl_tot = Cl_fada ;

%/ Cm_tot 
Cm_tot = Cm_fa + Cm_fade;

% Cn_tot
Cn_tot = Cn_fbdr + Cn_fada; 


% Total forces 
Xbar = qbar * Sref * CX_tot;
Ybar = qbar * Sref * CY_tot;
Zbar = qbar * Sref * CZ_tot;

% Total moments 
Lbar = Cl_tot * qbar * Sref * bref;
Mbar = Cm_tot * qbar * Sref * cref;
Nbar = Cn_tot * qbar * Sref * bref; 

F_Aerod_body = [(X(9)*Vbody(2) - X(8)*Vbody(3) + (Thrust - Xbar)/mass);(-X(9)*Vbody(1) + X(7)*Vbody(3) - Ybar/mass);(X(8)*Vbody(1) - X(7)*Vbody(2) - Zbar/mass)];

F_Aerod_ned = C_b2n*F_Aerod_body;

Fned = F_Aerod_ned - cross((2*w_ien + w_enn),[X(1);X(2);X(3)]) + g_ln;

% Forces along NED axis
DX(1) = Fned(1);

DX(2) = Fned(2);

DX(3) = Fned(3);

% Rate of change of euler angles w.r.t NED frame
DX(4) = X(7) + tan(X(5))*(X(8)*sin(X(4)) + X(9)*cos(X(4)));

DX(5) = X(8)*cos(X(4)) - X(9)*sin(X(4));

DX(6) = (X(8)*sin(X(4)) + X(9)*cos(X(4))) / cos(X(5));

% Angular acceleration
DX(7) = (C1*X(9) + C2*X(7))*X(8) + C3*Lbar + C4*(Nbar + X(8) * heng);

DX(8) = C5*X(7)*X(9) - C6*( X(7)^2 - X(9)^2) + C7*(Mbar - heng * X(9));

DX(9) = (C8*X(7) - C2*X(9))*X(8) + C4*Lbar + C9 * (Nbar + X(8) * heng);

% Ldot,ldot,hdot
DX(10) = X(1)/(R_N+X(12));

DX(11) = X(2)*sec(X(10))/(R_E+X(12));

DX(12) = -X(3);

AngV = [X(7) X(8) X(9)];
Acceleration_body = F_Aerod_body';
end