%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function evaluates the pseudo range measured by user        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pr_measured_L1, pr_measured_L2, pr_measured_L5] = eval_pr_measurement(gps_sat,sv_id,gps_time,true_user_pos_ecef, rcvr_clk_bias,true_sat_pos_ecef)

c = 2.99792458e8;

% [sv_x,sv_y,sv_z] = calc_sat_pos_ecef(gps_sat,gps_time,sv_id); % The true position of satellite based on ephemeris data and the gps time
% 
% true_sat_pos_ecef = struct('x',sv_x,'y',sv_y,'z',sv_z); % put this position in a struct

[sat_clk_drift,sat_clk_rel_error] = eval_sat_clock_offset(gps_sat,sv_id,gps_time); % satellite clock offset in seconds (clock drift + relativistic error)
        
time_str = gps_time + sat_clk_drift + sat_clk_rel_error; % The time of tranmission from satellite as embedded in transmitted code

D =  compute_distance(true_sat_pos_ecef,true_user_pos_ecef); % distance from user to satllite

user_pos_geodetic = ecef_to_latlong(true_user_pos_ecef); % user position in geodetic frame

[elev,azim,is_visible] = eval_el_az(user_pos_geodetic,true_user_pos_ecef,true_sat_pos_ecef); % gives elevation and azim in radians

slant_iono_delay_L1 = eval_delay_iono(user_pos_geodetic,elev,azim,gps_time); % calculate L1 ionospheric delay in seconds

slant_iono_delay_L2 = 1.6469*slant_iono_delay_L1; %calculate L2 ionospheric delay in seconds (from paper" a systems approach to ionospheric delay compenstaion")

slant_tropo_delay = eval_delay_tropo(elev,user_pos_geodetic.alt); % calculate tropospheric delay in seconds

rcvr_noise = randn; % Receiver error due to measurement noise

if rcvr_noise > 0.5 
    rcvr_noise = 0.5;
end;
if rcvr_noise <-0.5
    rcvr_noise = -0.5;
end;

time_rcvr_L1 = gps_time +  D/c + slant_tropo_delay + slant_iono_delay_L1 + rcvr_clk_bias; % Time of receiving as seen by user

time_rcvr_L2 = gps_time +  D/c + slant_tropo_delay + slant_iono_delay_L2 + rcvr_clk_bias; % Time of receiving as seen by user

pr_measured_L1 = c*(time_rcvr_L1 - time_str)+ rcvr_noise ; % computing the pseduo range which is measured by the receiver

pr_measured_L2 = c*(time_rcvr_L2 - time_str)+ rcvr_noise ; % computing the pseduo range which is measured by the receiver

A = (pr_measured_L1 - (2.546*pr_measured_L1 - 1.546*pr_measured_L2))*1575.42^2;

pr_measured_L5 = 2.546*pr_measured_L1 - 1.546*pr_measured_L2 + A/1176.45^2;


end