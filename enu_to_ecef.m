%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ENU to ECEF coordinate frame conversion function        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_ecef] = enu_to_ecef(ac_pos_enu, C_ne, ref_frame_pos_ecef)

position_ecef = [ref_frame_pos_ecef.x;ref_frame_pos_ecef.y;ref_frame_pos_ecef.z;] + C_ne*[ac_pos_enu.x; ac_pos_enu.y; ac_pos_enu.z];

x = position_ecef(1);
y = position_ecef(2);
z = position_ecef(3);

pos_ecef = struct('x',x,'y',y,'z',z);

end


