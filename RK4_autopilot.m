%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   RK4 INTEGRATION FUNCTION for autopilot   
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[Xnew] = RK4_autopilot(dt,X,Z,autopilotplant)

  XA = linspace(0,0,length(X));
  XB = XA;
  DX = acplant(X,Z,autopilotplant);
  XA = DX * dt;
  XB = X + 0.5 * XA;
  DX = acplant(X,Z,autopilotplant);
  Q =  DX * dt;
  XB = X + 0.5 * Q;
  XA = XA + 2.0 * Q;
  DX = acplant(X,Z,autopilotplant);
  Q = DX * dt;
  XB = X + Q;
  XA = XA + 2.0 * Q;
  DX = acplant(X,Z,autopilotplant);
  Xnew = X + (XA + DX * dt)/6.0;
end

function [DX] = acplant(X,Z,autopilotplant);
% H = eye(15,15); % measurement matrix
% K = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
Bw = [0.0210 -0.122;0.2090 0.5300;-0.0170 0.1640;zeros(12,2)];
[a,b,c,d] = ssdata(autopilotplant);
DX = (a*X' - b*[0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]*Z' + Bw*[15+randn;5+randn])';

end