%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Distance Computation Function        %
%   Author: Saurav Agarwal   %
%   Date: January 1, 2011  %
%   Dept. of Aerospace Engg., IIT Bombay, Mumbai, India %
%   All i/o units are specified in brackets %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:     
%        1. distance: (m) 
% Inputs:
%        1. pos1/pos2: pos of 2 points in ECEF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [distance] = compute_distance(pos1,pos2);

distance =  sqrt((pos1.x-pos2.x)^2 + (pos1.y-pos2.y)^2+(pos1.z-pos2.z)^2 );

end