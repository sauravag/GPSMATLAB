%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Atmosphere Function   %
%   Author: Saurav Agarwal   %
%   Date: January 1, 2011  %
%   Dept. of Aerospace Engg., IIT Bombay, Mumbai, India %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FUNCTION TO CALCULATE 'qbar'
function [qbar,rho] = atmosphere(alt,Vt)

    %** CONSTANT PARAMETERS **
g = 9.86; % [m/s^2]
temp_lapserate = 0.005; % [K/m]
temp_sealevel  = 300; % Temperature at sealevel[K]
pressure_sealevel = 101325; % Pressure at sea level [pa]

temp_alt = temp_sealevel - (alt * temp_lapserate);
pressure_alt = pressure_sealevel * ((temp_sealevel/temp_alt)^(-g/(temp_lapserate*287.0)));
rho = pressure_alt / (287.0*temp_alt);
qbar = 0.5 * rho * Vt * Vt;
end