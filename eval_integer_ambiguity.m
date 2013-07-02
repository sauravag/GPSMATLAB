%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function evaluates the integer ambiguities        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IntegerAmbiguities,AA,clkbias_u,clkbias_r] = eval_integer_ambiguity(no_of_measurement_epochs,time_of_epochs,countme,gps_sat,visible_sats_id,ref_station_ecef,ref_station,initial_user_pos_estimate_storage,true_user_pos_geodetic_storage,true_user_pos_ecef_storage,ibeacon1_geo,ibeacon1_ecef,ibeacon2_geo,ibeacon2_ecef,initial_clk_bias_u,initial_clk_bias_r);

for i=1:countme

    gps_time = time_of_epochs(i);

    [deltaPhi,Sk,SmoothPhase,clkbias_u,clkbias_r] = collect_phase_data(gps_sat,gps_time,visible_sats_id,ref_station_ecef,ref_station,true_user_pos_geodetic_storage(i),true_user_pos_ecef_storage(i), initial_user_pos_estimate_storage(i), ibeacon1_geo,ibeacon1_ecef,ibeacon2_geo,ibeacon2_ecef,initial_clk_bias_u,initial_clk_bias_r);

    AA(:,i) = SmoothPhase; % Each column contains the smoothed phase (1st row contains N1 and so on)

end;
[rows,columns] = size(AA)
for j = 1:rows
    IntegerAmbiguities(j) = (1/0.19)*sum(AA(j,:))/columns ; 
end;    
        
end