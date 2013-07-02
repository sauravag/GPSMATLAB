%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Ionospheric Delay Calculation
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Function for   computing an Ionospheric range correction for the *
%      GPS L1 frequency from the parameters broadcasted in the GPS      *
%      Navigation Message.                                              *
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      References:                                                      *
%      Klobuchar, J.A., (1996) "Ionosphercic Effects on GPS", in        *
%        Parkinson, Spilker (ed), "Global Positioning System Theory and *
%        Applications, pp.513-514.                                      *
%      ICD-GPS-200, Rev. C, (1997), pp. 125-128                         *
%      NATO, (1991), "Technical Characteristics of the NAVSTAR GPS",    *
%        pp. A-6-31   -   A-6-33                                        *
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input :                                                            *
%        user_pos_geodetic
%        El : elevation in radians                                       %                                           
%        A : azimuth in radians                                         %
%        GPS_Time      : Time of Week                           (sec)   *
%        Alfa(4)       : The coefficients of a cubic equation           *
%                        representing the amplitude of the vertical     *
%                        dalay (4 coefficients - 8 bits each)           *
%        Beta(4)       : The coefficients of a cubic equation           *
%                        representing the period of the model           *
%                        (4 coefficients - 8 bits each)                 *
%    Output:                                                            *
%       Iono_delay        : Ionospheric slant range correction in sec   *
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Iono_delay]= eval_delay_iono(user_pos_geodetic,El,A,gps_time);

    % Alpha and Beta parameters for Iono model taken from
    % ftp://cddis.gsfc.nasa.gov/pub/gps/data/daily
    Alpha =  [0.1118E-07 0.0000E+00 -0.5960E-07 0.0000E+00];
    Beta  =  [0.9011E+05 0.0000E+00  -0.1966E+06   0.0000E+00];


    lat_user = user_pos_geodetic.lat/3.14159; % Convert to semi-circles
    
    long_user = user_pos_geodetic.long/3.14159; % Convert to semi-circles
     
    El = El/3.14159; % Convert to semi-circles
    
    %Calculate the earth centered angle
    Psi = 0.0137/(El+0.11)-0.022; %in semi-circles
    
    %Compute the subionospheric latitude
    Phi_subiono = lat_user + Psi*cos(A);% semi-circles;
    
    if Phi_subiono > 0.416
        
        Phi_subiono = 0.416;
    end;
    if Phi_subiono < -0.416
        Phi_subiono = -0.416;
    end;
    
    %Compute the subionospheric longitude
    lambda_subiono = long_user + Psi*sin(A)/cos(Phi_subiono);
    
    %compute the geomagnetic latitude
    Phi_m = Phi_subiono + 0.064*cos(lambda_subiono - 1.617);
    
    % find the local time
    t = 4.32*10000*lambda_subiono + gps_time;
    
    if t > 86400
        t = t - 86400;
    end;
    if t<0
        t = t + 86400;
    end;
    
    %compute the slant delay factor
    F = 1.0 + 16.0*(0.53-El)^3;
    
    Beta_factor=0;
    Alpha_factor = 0;
    for n=1:4
        Beta_factor =  Beta_factor + Beta(n)*Phi_m^(n-1);
        Alpha_factor = Alpha_factor + Alpha(n)*Phi_m^(n-1);
    end;  
    
    if Beta_factor < 72000                                    
        Beta_factor=72000;
    end;
    
    if Alpha_factor < 0 
        Alpha_factor = 0;
    end
    % Evaluate the delay 
    
    x = 2*3.14159*(t - 50400)/(Beta_factor);
    
    if abs(x) > 1.57
        
        Iono_delay = F*5e-9;
    end;
    if abs(x) <= 1.57
        
        Iono_delay = F*(5e-9+Alpha_factor*(1-(x^2)/2+(x^4)/24));
    end;
    
end

