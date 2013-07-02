function [pos_ned]  = eval_ned_to_ecef(ref_station,  ref_station_ecef, user_pos_ned)
lat = ref_station.lat;
long = ref_station.long;

xu = user_pos_ned.x;
yu = user_pos_ned.y;
zu = user_pos_ned.z;

xr = ref_station_ecef.x;
yr = ref_station_ecef.y;
zr = ref_station_ecef.z;

C_ECEF2NED = [-sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat);-sin(long) cos(long) 0; -cos(lat)*cos(long) -cos(lat)*sin(long) -sin(lat)];

r_ecef = C_ECEF2NED'*[xu;yu;zu];

pos_ned = struct('x', r_ecef(1)+xr,'y',r_ecef(2)+yr,'z',r_ecef(3)+zr);
end
