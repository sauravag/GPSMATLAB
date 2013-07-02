%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Clock bias function oscillator type: TCXO    
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rcvr_clock_error] = rcvr_clk_model(initial_bias,randomfactor)

  
 rcvr_clock_error = initial_bias + 2.5e-6 + 1e-6/(2*32768*365*24*3600) + 5e-8*randomfactor ;
 
end