%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global Positioning System Simulation Matlab Tool        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function contains the location of the indian reference stations

function [list_inres] = inres_pos_data()

inres_1 = struct('lat',28.55*dtr,'long',77.1*dtr,'alt',216); %delhi
inres_2 = struct('lat',32.683*dtr,'long',74.233,'alt',314); %jammu
inres_3 = struct('lat',23.066*dtr,'long',38.066*dtr,'alt',58); %ahmedabad
inres_4 = struct('lat',26.25*dtr,'long',73.02*dtr,'alt',219); %jodhpur
inres_5 = struct('lat',26.75*dtr,'long',77.33*dtr,'alt',125); % lucknow
inres_6 = struct('lat',23.28*dtr,'long',77.33*dtr,'alt',524);%bhopal
inres_7 = struct('lat',26.67*dtr,'long',88.33*dtr,'alt',126); %bagdogra
inres_8 = struct('lat',26.1*dtr,'long',91.5833*dtr,'alt',49);%guwahai
inres_9 = struct('lat',23.833*dtr,'long',92.62*dtr,'alt',405); % aizwal
inres_10 = struct('lat',27.5*dtr,'long',95*dtr,'alt',110);%dibrugarh
inres_11 = struct('lat',21.17*dtr,'long',81.73*dtr,'alt',317);%raipur
inres_12 = struct('lat',22.65*dtr,'long',88.433*dtr,'alt',5); %kolkata
inres_13 = struct('lat',19.083*dtr,'long',72.86*dtr,'alt',11);%mumbai
inres_14 = struct('lat',17.28*dtr,'long',78.42*dtr,'alt',617);%hyderabad
inres_15 = struct('lat',17.72*dtr,'long',83.22*dtr,'alt',5); %vishakhapatnam
inres_16 = struct('lat',13.20*dtr,'long',77.7*dtr,'alt',915); %bangalore
inres_17 = struct('lat',13*dtr,'long',9.15*dtr,'alt',16); %chennai
inres_18 = struct('lat',10.82*dtr,'long',72.17*dtr,'alt',4); %agatti 
inres_19 = struct('lat',8.5*dtr,'long',76.92*dtr,'alt',4); %trivandrum
inres_20 = struct('lat',11.633*dtr,'long',92.72*dtr,'alt',4);%port blair

list_inres = [inres_1,inres_2,inres_3,inres_4inres_5,inres_6,inres_7,inres_8,inres_9,inres_10,inres_11,inres_12,inres_13,inres_14,inres_15,inres_16,inres_17,inres_18,inres_19,inres_20];
end













