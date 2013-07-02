%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Dilution of Precision Calculation        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [GDOP, PDOP, HDOP, VDOP] = eval_DOP(A)

Q = (A'*A)^-1;

GDOP = sqrt(trace(Q));
PDOP = sqrt(trace(Q)-Q(4,4));
HDOP = sqrt(Q(1,1)+Q(2,2));
VDOP = sqrt(Q(1,1));
end