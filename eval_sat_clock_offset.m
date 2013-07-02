%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate Satellite Offest Clock Error        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:"GPS Theory and application",edited by B.Parkinson,J.Spilker, 
%Ofset Model:                                                             
%dTclk_Ofset=af0+af1*(T-Toc)+af2*(T-Toc)^2+.....                          
%af :(1/Sec^i)    =>Matrix of Coeeficient for satellite offset            
%Ttr:(Sec)        => Time of transmission                                 
%Toc:(Sec)        => Sv Clock refernce time                               
%dTclk_Ofset:(Sec)=> Sv Clock offset time                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sat_clk_drift,sat_clk_rel_error]= eval_sat_clock_offset(gps_sat,sv_id,t_sv)

    mu = 3.986005e14;
    c = 2.99792458e8;	
    F = -2*sqrt(mu)/c^2;
    delta_n = 4.908419e-9;
    t_oe = 147456.0000;
    
    n_0 = sqrt(mu/(gps_sat(sv_id).sqrt_a)^6);		% (rad/s)
    t_k=t_sv-t_oe;									% Time from eph ref epoch (s)
    n = n_0 + delta_n;                              % Corrected mean motion (rad/s)
    M_k = gps_sat(sv_id).M_0+n*t_k;					% Mean anomaly (rad/s)

    %  Perform Newton-Raphson solution for E_k estimate
    E_k= newton_raphson(gps_sat(sv_id).e,M_k);		% Eccentric anomaly ESTIMATE for computing delta_tr
    
    sat_clk_rel_error = F*gps_sat(sv_id).e*gps_sat(sv_id).sqrt_a*sin(E_k); % relativistic clock correction in seconds
    
    sat_clk_drift = gps_sat(sv_id).af0 + gps_sat(sv_id).af1*(t_sv - t_oe); %drift in in seconds
    
end

    