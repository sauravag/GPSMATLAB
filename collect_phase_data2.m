%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      GPS Signal Phase Measurement Function        %
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
%            1. "Global Positioning System: Theory and Applications Vol II, Bradford W. Parksinson, Pg 435-437" 
% Outputs:     
%        1. deltaPhi: single difference phase measurements between user and ref 
%        2. Scapk: geometry matrix
%        3. SmoothPhase: smoothed pseudorange
%        4. clbias_u: user clock bias (m) 
%        5. clkbias_r: reference clock bias (m)
% Inputs:
%        1. gps_sat: array containing ephemeris data of gps satellites
%        2. gps_time: gps time (s)
%        3. sv_id: id number of gps satellite for which to calculate p/v
%        4. visible_sats_id:
%        5. ref_station_ecef
%        6. ref_station: reference station pos in 
%        7. true_user_pos_geodetic: true user pos in 
%        8. true_user_pos_ecef: true ECEF pos (m)
%        9. initial_user_pos_estimate: initial pos estimate in ECEF (m)
%       10. ibeacon1_geo/ibeacon1_ecef: integrity beacon 1 position in geodetic and ECEF (m) coordinates respctly
%       11. ibeacon2_geo/ibeacon2_ecef: integrity beacon 2 position in geodetic and ECEF (m) coordinates respctly
%       12. initial_bias_u/initial_bias_r: initial clock bias of user and reference reciever respctly (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [deltaPhi,Scapk,SmoothPhase] = collect_phase_data2(gps_sat,gps_time,visible_sats_id,ref_station_ecef,ref_station,true_user_pos_geodetic,true_user_pos_ecef,initial_user_pos_estimate,ibeacon1_geo,ibeacon1_ecef,ibeacon2_geo,ibeacon2_ecef,clkbias_u,clkbias_r)

c = 2.99792458e8;% speed of light (m/s)
f_L1 = 1575.42e6; %frequency of L1 carrier wave Hz
xu = initial_user_pos_estimate.x;
yu = initial_user_pos_estimate.y;
zu = initial_user_pos_estimate.z;

xutrue = true_user_pos_ecef.x;
yutrue = true_user_pos_ecef.y;
zutrue = true_user_pos_ecef.z;

user_ecef_pos_estimate = struct('x',xu,'y',yu,'z',zu);

% coordinates of reference station & integrity beacons
xr = ref_station_ecef.x;
yr = ref_station_ecef.y;
zr = ref_station_ecef.z;

xib1 = ibeacon1_ecef.x;
yib1 = ibeacon1_ecef.y;
zib1 = ibeacon1_ecef.z;

xib2 = ibeacon2_ecef.x;
yib2 = ibeacon2_ecef.y;
zib2 = ibeacon2_ecef.z;

dim = length(visible_sats_id);

lambda_L1 = 0.19; % metres

range_ru_true = compute_distance(true_user_pos_ecef,ref_station_ecef); % the true range between the reference station and satellite
range_uib1_true = compute_distance(true_user_pos_ecef,ibeacon1_ecef);
range_uib2_true = compute_distance(true_user_pos_ecef, ibeacon2_ecef);
range_rib1_true = compute_distance(ref_station_ecef,ibeacon1_ecef);
range_rib2_true = compute_distance(ref_station_ecef, ibeacon2_ecef);

x = [xu-xr;yu-yr;zu-zr]; % relative vector from ref to user

randomfactor_u = randn; % gaussian noise factor in user clock
randomfactor_r = randn; % gaussian noise factor in reference clock
randomfactor_ib1 = randn; % gaussian noise factor in user clock
randomfactor_ib2 = randn; % gaussian noise factor in reference clock
  
for k = 1:dim
    
     sv_id = visible_sats_id(k);
            
     [xs(k),ys(k),zs(k)] = calc_sat_pos_ecef(gps_sat,gps_time,sv_id); % The true position of satellite based on ephemeris data and the time embedded in message

     true_sat_pos_ecef = struct('x',xs(k),'y',ys(k),'z',zs(k)); 
     
     r_rs(k) = compute_distance(true_sat_pos_ecef,ref_station_ecef); % the true range between the reference station and satellite
     
     r_us(k) = compute_distance(true_sat_pos_ecef,true_user_pos_ecef); 
    
     N = k;
    if N >3
         N = k - 3;
     end;
     if N > 6
         N = k - 3;
     end;
     if N==1
         N = 0;
     end;
     
     trueambiguity(k) = (-1)^k*N*1057;
     
     s = [(xs(k)-xr);(ys(k)-yr);(zs(k)-zr)]/r_rs(k);%[(xs(k)-(xu+xr)/2);(ys(k)-(yu+yr)/2);(zs(k)-(zu+zr)/2)]/r_rs(k) ; % the line of sight vector from ref station to satelite k
     
     psi_ru(k) = r_us(k) - r_rs(k) + (-1)^k*N*1057*lambda_L1 + (clkbias_u -clkbias_r)*c + randn/100; % the single difference between the phase at reference station and user for satellite k
     
     [pr_measured_L1_u, pr_measured_L2, pr_measured_L5] = eval_pr_measurement(gps_sat,sv_id,gps_time,true_user_pos_ecef,clkbias_u,true_sat_pos_ecef);
     
     [pr_measured_L1_r, pr_measured_L2, pr_measured_L5] = eval_pr_measurement(gps_sat,sv_id,gps_time,ref_station_ecef,clkbias_r,true_sat_pos_ecef);
     
     rho_ru(k) = pr_measured_L1_u - pr_measured_L1_r ; % single difference code measurement 
     
     delta_psi_ru(k) = psi_ru(k)+ s'*x;
     
     Scap(k,:) = [-s' 1];
end;

e1 = [xib1-xu;yib1-yu;zib1-zu]/norm([xib1-xu;yib1-yu;zib1-zu]);
e2 = [xib2-xu;yib2-yu;zib2-zu]/norm([xib2-xu;yib2-yu;zib2-zu]);

e1dash = [-e1' 1];
e2dash = [-e2' 1];

Scapk = [Scap;e1dash;e2dash];

p1 = [xib1-xr;yib1-yr;zib1-zr];% rel vector from ref to ib1
p2 = [xib2-xr;yib2-yr;zib2-zr];% rel vector from ref to ib2

phi_ruib1 = (range_uib1_true-range_rib1_true) + 3100*lambda_L1 + (clkbias_u -clkbias_r)*c + randn/100; %the single difference between the phase at user and ref station for integrity beacon 1 

rcvr_noise = randn/2; % Receiver error due to measurement noise

if rcvr_noise > 0.5 
    rcvr_noise = 0.5;
end;
if rcvr_noise <-0.5
    rcvr_noise = -0.5;
end;

rho_ruib1 = (range_uib1_true-range_rib1_true) + (clkbias_u -clkbias_r)*c + rcvr_noise ;

phi_ruib2 = (range_uib2_true-range_rib2_true) + 2200*lambda_L1 + (clkbias_u -clkbias_r)*c  + randn/100; %the single difference between the phase at user and ref station for integrity beacon 2
rcvr_noise = randn/2; % Receiver error due to measurement noise

if rcvr_noise > 0.5 
    rcvr_noise = 0.5;
end;
if rcvr_noise <-0.5
    rcvr_noise = -0.5;
end;
rho_ruib2 = (range_uib2_true-range_rib2_true) + (clkbias_u -clkbias_r)*c  + rcvr_noise;

delta_phi_ruib1 = phi_ruib1 - (norm(p1-x) - norm(p1));
delta_phi_ruib2 = phi_ruib2 - (norm(p2-x) - norm(p2));

deltaPhi = [delta_psi_ru';delta_phi_ruib1;delta_phi_ruib2];
SDPhi = [psi_ru';phi_ruib1;phi_ruib2];
SDRho = [rho_ru';rho_ruib1;rho_ruib2];
SmoothPhase =   SDPhi - SDRho;
actualambiguity = [trueambiguity';3100;2200];
end   