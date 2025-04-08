function [d_x_1, d_x_2, d_x_3, d_x_4] = dip_dynamics(x_1, x_2, x_3, x_4, u)
  % dynamics of the dip from
  % https://link.springer.com/content/pdf/10.1134/S1064230706030014.pdf
  % with the state 
  %   x_1 = theta_1
  %   x_2 = theta_2
  %   x_3 = d_theta_1
  %   x_4 = d_theta_2
  % parameters
  a_11 = 1;
  a_12 = 1;
  a_22 = 1;
  b_1 = 1;
  b_2 = 1;
  A = [a_11 a_12*cos(x_2-x_1); 
       a_12*cos(x_2-x_1) a_22];
  F = [0 -a_12*sin(x_2-x_1);
       a_12*sin(x_2-x_1)];
  B = [-b_1 0; 0 -b_2]
  % change c depending on analysis
  c = [1;0];
  % c = [-1;1];
  % compute derivative dynamics
  dd_theta = inv(A)*(c*u-F*[x_3^2;x_4^2]-B*[sin(x_1);sin(x_2)])
  d_x_1 = x_3;
  d_x_2 = x_4;
  d_x_3 = dd_theta(1);
  d_x_4 = dd_theta(2);
end
