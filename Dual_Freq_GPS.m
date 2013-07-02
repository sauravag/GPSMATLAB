%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Dual frequency GPS Reciever Emulation Function        %
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   All i/o units are specified in brackets %
%   Conventions:
%                1. WGS-84 system used for geodetic model
%                2. Geodetic Coordinates are used in the form lat(rad)/long(rad)/alt(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
%            1. Global Positioning System: Theory and Applications Vol I,Bradford W. Parksinson, 
% Outputs:     
%        1. user_pos_gps: reciver pos as calculated from gps measurements ECEF coordinates (m)
%        2. optimum_sv_ids: list of satellites with lowest GDOP 
%        3. DOP: array containing GDOP,PDOP,HDOP,VDOP
%        4. Velocity_Ecef: Velocity of reciever in ECEF frame (m)
%        5. rcvr_clk_bias: clock bias (m)
% Inputs:
% true_user_pos_ecef,initial_user_pos_estimate,inres_pos_data,estimate_user_vel_ecef,true_user_vel_ecef
%        1. initial_bias: initial clock bias (s)
%        2. GPSMODE: all in view/4 satellite mode
%        3. gps_sat: gps satellite ephemerides
%        4. gps_time: (s)
%        5. visible_sats_id:list of visible gps satellites 
%        6. true_user_pos_ecef: true ECEF pos (m)
%        7. initial_user_pos_estimate: initial pos estimate in ECEF (m)
%        8. inres_pos_data: indian GAGAN reference station position data
%        9. estimate_user_vel_ecef: estimtae of user velocity (m/s) ECEF frame
%       10. true_user_vel_ecef: true user velocity (m/s) ECEF frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [user_pos_gps,optimum_sv_ids,DOP,Velocity_Ecef,rcvr_clk_bias] = Dual_Freq_GPS(initial_bias,GPSMODE,gps_sat,gps_time,visible_sats_id,true_user_pos_ecef,initial_user_pos_estimate,inres_pos_data,estimate_user_vel_ecef,true_user_vel_ecef)

c = 2.99792458e8;% speed of light (m/s)

xu = initial_user_pos_estimate.x;
yu = initial_user_pos_estimate.y;
zu = initial_user_pos_estimate.z;

user_ecef_pos_estimate = struct('x',xu,'y',yu,'z',zu);

N = length(visible_sats_id);

[optimum_sv_ids, GDOP,PDOP,HDOP,VDOP] = select_optimum_sats(gps_sat,gps_time,visible_sats_id,true_user_pos_ecef,initial_user_pos_estimate);
DOP = [GDOP PDOP HDOP VDOP];
delta_x = 10;
delta_y = 10;
delta_z = 10;
Cb = 0;
rcvr_clk_bias = rcvr_clk_model(initial_bias,randn);
%sbas_correction = 0;%eval_sbas_correction(gps_time,gps_sat,inres_pos_data,true_user_pos_ecef); % correction in seconds provided by aigmentation system
% Initialise arrays for storage
% sat_clk_drift = zeros(1,4);
% sat_clk_rel_error =zeros(1,4);
% xs = zeros(1,4);
% ys = zeros(1,4);
% zs = zeros(1,4);
% pr_measured_L1 = zeros(1,4);
% pr_measured_L2 = zeros(1,4);
% pr_measured_L5 = zeros(1,4);
% Vsat_ECEF = zeros(4,3);
% pr_ionofree = zeros(1,4);
% Gtrue =  zeros(4,3);
%  L = zeros(1,4);
%  A = zeros(4,4);
%  AA = zeros(4,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All in view mode %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 if GPSMODE == 1
    
     for k = 1:N
     
     sv_id = visible_sats_id(k);%optimum_sv_ids(k);
     
     [sat_clk_drift(k),sat_clk_rel_error(k)] = eval_sat_clock_offset(gps_sat,sv_id,gps_time); % satellite clock offset in seconds (clock drift + relativistic error)
     
     [xs(k),ys(k),zs(k),Vsat_ECEF(k,:)] = calc_sat_pos_ecef(gps_sat,gps_time,sv_id); % The position of satellite based on ephemeris data and the time embedded in message

     computed_sat_pos_ecef(k) = struct('x',xs(k),'y',ys(k),'z',zs(k));
         
     [pr_measured_L1(k), pr_measured_L2(k), pr_measured_L5] = eval_pr_measurement(gps_sat,sv_id,gps_time,true_user_pos_ecef, rcvr_clk_bias,computed_sat_pos_ecef(k)); 
     
     pr_ionofree(k) = pr_measured_L1(k)*2.546 - pr_measured_L2(k)*1.546 + c*(sat_clk_drift(k)+sat_clk_rel_error(k)); % make it ionofree from page 166 mishra and enge and correct it for satellite clock error
         
     dtrue =  compute_distance(computed_sat_pos_ecef(k),true_user_pos_ecef);
     
     Gtrue(k,:) = [(-(xs(k) - true_user_pos_ecef.x)/dtrue) (-(ys(k) - true_user_pos_ecef.y)/dtrue) (-(zs(k) - true_user_pos_ecef.z)/dtrue)];
     end; 
     
    while abs(delta_x)> 1e-3 && abs(delta_y)> 1e-3 && abs(delta_z)> 1e-3

     for k = 1:N
        
         d = compute_distance(computed_sat_pos_ecef(k),user_ecef_pos_estimate);
         
         pr_estimated =  d + Cb ; % what the receiver thinks is true pseduo range

         L(k) = pr_ionofree(k) - pr_estimated ; %delta rho
         
         A(k,:) = [(-(xs(k) - xu)/d) (-(ys(k) - yu)/d) (-(zs(k) - zu)/d) 1];

     end; 

     error = (A'*A)^-1*A'*L';

     Cb = Cb + error(4);

     delta_x = error(1);

     delta_y = error(2);

     delta_z = error(3);

     xu = xu + delta_x ;

     yu = yu+ delta_y ;

     zu = zu + delta_z ;

     user_ecef_pos_estimate = struct('x',xu,'y',yu,'z',zu);
    end
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Optimal satellite geometry mode %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if GPSMODE == 2

      for k = 1:4
     
          sv_id = optimum_sv_ids(k);%optimum_sv_ids(k);

         [sat_clk_drift(k),sat_clk_rel_error(k)] = eval_sat_clock_offset(gps_sat,sv_id,gps_time); % satellite clock offset in seconds (clock drift + relativistic error)

         [xs(k),ys(k),zs(k),Vsat_ECEF(k,:)] = calc_sat_pos_ecef(gps_sat,gps_time,sv_id); % The position of satellite based on ephemeris data and the time embedded in message

         computed_sat_pos_ecef(k) = struct('x',xs(k),'y',ys(k),'z',zs(k));

         [pr_measured_L1(k), pr_measured_L2(k), pr_measured_L5] = eval_pr_measurement(gps_sat,sv_id,gps_time,true_user_pos_ecef, rcvr_clk_bias,computed_sat_pos_ecef(k)); 

         pr_ionofree(k) = pr_measured_L1(k)*2.546 - pr_measured_L2(k)*1.546 + c*(sat_clk_drift(k)+sat_clk_rel_error(k)); % make it ionofree (from page 166 mishra and enge) and correct for satellite clock error

         dtrue =  compute_distance(computed_sat_pos_ecef(k),true_user_pos_ecef);
         
         Gtrue(k,:) = [(-(xs(k) - true_user_pos_ecef.x)/dtrue) (-(ys(k) - true_user_pos_ecef.y)/dtrue) (-(zs(k) - true_user_pos_ecef.z)/dtrue)];
     end; 
     
    while abs(delta_x)> 1e-3 && abs(delta_y)> 1e-3 && abs(delta_z)> 1e-3

     for k = 1:4

         d = compute_distance(computed_sat_pos_ecef(k),user_ecef_pos_estimate); % what the receiver thinks is true pseduo range

         pr_estimated = d + Cb;
         
         L(k) = pr_ionofree(k) - pr_estimated ; % delta rho
         
         A(k,:) = [(-(xs(k) - xu)/d) (-(ys(k) - yu)/d) (-(zs(k) - zu)/d) 1];

         
     end; 

     error = A^-1*L';

     Cb = Cb + error(4);

     delta_x = error(1);

     delta_y = error(2);

     delta_z = error(3);

     xu = xu + delta_x ;

     yu = yu+ delta_y ;

     zu = zu + delta_z ;

     user_ecef_pos_estimate = struct('x',xu,'y',yu,'z',zu);
    end
 end
user_pos_gps = user_ecef_pos_estimate;

Velocity_Ecef = calc_user_velocity_ecef(A,Gtrue,Vsat_ECEF,estimate_user_vel_ecef,true_user_vel_ecef);

end

% Velocity estimation using doppler shift
% Page 411 Parkinson Vol.I
function [Velocity_ECEF] = calc_user_velocity_ecef(A,Gtrue,Vsat_ECEF,estimate_user_vel_ecef,true_user_vel_ecef);

c = 2.99792458e8;% speed of light (m/s)
delta_vx = 10;
delta_vy = 10;
delta_vz = 10;
f = c*1e-6/(365*24*3600); % clock bias in m/s
N = length(A);
noise = randn;
clkdrift_estimate = 0;
while abs(delta_vx)> 1e-3 && abs(delta_vy)> 1e-3 && abs(delta_vz)> 1e-3
    
    for i = 1:N

        rho_dot(i) =  dot((Vsat_ECEF(i)'-true_user_vel_ecef'),Gtrue(i,1:3)') + f + noise; % Measured pseudorange rate

        rho_cap_dot(i) = dot((Vsat_ECEF(i)'- estimate_user_vel_ecef'),A(i,1:3)') + clkdrift_estimate ; % estimated pseudo range rate

        delta_rho_dot(i) = rho_cap_dot(i) - rho_dot(i); 
        
        
    end
    
    error =  (A'*A)^-1*A'*delta_rho_dot'; % least squares
    
    clkdrift_estimate = clkdrift_estimate + error(4);
    
    delta_vx = error(1);
 
    delta_vy = error(2);

    delta_vz = error(3);

    estimate_user_vel_ecef = estimate_user_vel_ecef + [error(1) error(2) error(3)];

end;

    Velocity_ECEF = estimate_user_vel_ecef; % Estimated user velocity in ECEF frame
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    