function [d_x_1, d_x_2, d_x_3, d_x_4] = dip_dynamics(x_1, x_2, x_3, x_4)
  % dynamics of the dip from
  % https://link.springer.com/content/pdf/10.1134/S1064230706030014.pdf
  % with the state 
  %   x_1 = theta_1
  %   x_2 = d_theta_1
  %   x_3 = theta_2
  %   x_4 = d_theta_2
  % TODO: need to rearrange the dynamics in the paper to get d_x_2 and d_x_4
  d_x_1 = x_2;
  d_x_2 = 0;
  d_x_3 = x_4;
  d_x_4 = 0;
end
