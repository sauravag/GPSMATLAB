% FUNCTION TO CALCULATE THRUST
function [Thrust] = eval_Thrust(rho,Vt,dth);
 Pi =  3.141592654;
  J = Vt / (dth * 0.3);
  Thrust = 33 * (4/(Pi^2)) * rho * (dth ^ 2) * (0.3^4) * (-0.0948 * (J^2) + 0.058 * J + 0.0761 );
end