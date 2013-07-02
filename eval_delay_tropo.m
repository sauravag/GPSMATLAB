%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Tropospheric Delay Calculation
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                           *
% Reference:
%  "Global Positioning System, Mishra & Enge", pg 172            *
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%Input
%        T_amb:Degree Celsius =>At reciever antenna location
%        P_amb:kPa =>At reciever antenna location
%        P_vap:kPa =>Water vapore pressure at reciever antenna location
%        E => Elevation angle in radians from user to satellite 
%Output
%        slant tropo_delay in seconds

function [tropo_delay] = eval_delay_tropo(El,alt);

    [T_o,P_o] = atmosphere(alt); % local conditions
    P_vap = 20; % local vapour pressure  (milli bar)
    c = 2.99792458e8;% speed of light (m/s)
    hd = 43000; 
    hw = 12000;

    T_d = 77.6e-6*(P_o/T_o)*(hd/5);
    
    T_w = 0.373*(P_vap/T_o^2)*(hw/5);
    
    m_d = 1/(sin(El) + 0.00143/(tan(El)+0.0445));
    
    m_w = 1/(sin(El)+ 0.00035/(tan(El)+0.017));
    
    %Troposhpheric Delay 
    tropo_delay = (T_d*m_d + T_w*m_w)/c; % in seconds
end

function [temp, pres] = atmosphere(alt)

% CONSTANT PARAMETERS 
g = 9.86; % [m/s^2]
temp_lapserate = 0.005; % [K/m]
temp_sealevel  = 308; % Temperature at sealevel[K]
pressure_sealevel = 101325; % Pressure at sea level [pa]
temp_alt = temp_sealevel - (alt * temp_lapserate);
pressure_alt = pressure_sealevel * ((temp_sealevel/temp_alt)^(-g/(temp_lapserate*287.0)));
temp = temp_alt;
pres = pressure_alt*10^-2; % pressure in milli bar

end
    