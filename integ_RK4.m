%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   RK4 INTEGRATION FUNCTION        
%   Author: Saurav Agarwal   
%   Email:  saurav6@gmail.com
%   Date:   January 1, 2011  
%   Place:  Dept. of Aerospace Engg., IIT Bombay, Mumbai, India 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Xnew,DX,Acceleration_body,AngV] = integ_RK4(dt,X,dth,de,da,dr);

  [DX,Acceleration_body(1,:),AngV(1,:)] = eval_state_derivative(X,dth,de,da,dr);
  XA = DX * dt;
  XB = X + 0.5 * XA;
  [DX,Acceleration_body(2,:),AngV(2,:)] = eval_state_derivative(XB,dth,de,da,dr);
  Q =  DX * dt;
  XB = X + 0.5 * Q;
  XA = XA + 2.0 * Q;
  [DX,Acceleration_body(3,:),AngV(3,:)] = eval_state_derivative(XB,dth,de,da,dr);
  Q = DX * dt;
  XB = X + Q;
  XA = XA + 2.0 * Q;
  [DX,Acceleration_body(4,:),AngV(4,:)] = eval_state_derivative(XB,dth,de,da,dr);
  Xnew = X + (XA + DX * dt)/6.0;
  
end