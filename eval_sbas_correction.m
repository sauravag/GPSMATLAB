function [sbas_correction] = eval_sbas_correction(gps_time,gps_sat,inres_pos_data,user_pos_ecef)

c = 2.99792458e8;% speed of light (m/s)

for i = 1:20
    
    inres_ecef = latlong_to_ecef(inres_pos_data(i));
    d(i) = compute_distance(user_pos_ecef,inres_ecef);
end;

temp = min(d);

for j = 1:20
    
    if d(j) == min(d)
        closest_inres = j;
    end;
end;

closest_inres_pos_ecef = latlong_to_ecef(inres_pos_data(closest_inres));
closest_inres_pos_geodetic = inres_pos_data(closest_inres);

N = 0; %no of visible sats

for sv_id=1:31 % space vehicle number
            
        [sv_x,sv_y,sv_z] = calc_sat_pos_ecef(gps_sat,gps_time,sv_id); % The true position of satellite based on ephemeris data and the gps time

        true_sat_pos_ecef = struct('x',sv_x,'y',sv_y,'z',sv_z);
         
        [elev,azim,is_visible] = eval_el_az(closest_inres_pos_geodetic,closest_inres_pos_ecef,true_sat_pos_ecef); % gives elevation and azim in radians
             
        if is_visible == 1
            
            N = N +1;
            iono_delay(N) = eval_delay_iono(closest_inres_pos_geodetic,elev,azim,gps_time);
            tropo_delay(N) = eval_delay_tropo(elev,closest_inres_pos_geodetic.alt);
                                          
        end;   
        
end;

iono_correction = mean(iono_delay);% + std(iono_delay);
tropo_correction = mean(tropo_delay);% + std(tropo_delay);
sbas_correction = c*(iono_correction + tropo_correction);% delay correction in metres

end
        