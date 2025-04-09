function dx = dip_dynamics(t,x)
  % dynamics of the dip from
  % https://link.springer.com/content/pdf/10.1134/S1064230706030014.pdf
  % with the state 
  %   x_1 = theta_1
  %   x_2 = theta_2
  %   x_3 = d_theta_1
  %   x_4 = d_theta_2
  % controller
  u = 0;
  % parameters
  m_1 = 1;
  m_2 = 1;
  I_1 = 1;
  I_2 = 1;
  l = 2;
  r_1=1;
  r_2=1;
  g=9.81;
  a_11 = I_1+m_2*l^2;
  a_12 = m_2*r_2*l;
  a_22 = I_2;
  b_1 = (m_1*r_1+m_2*l)*g;
  b_2 = (m_2*r_2*g);
  A = [a_11 a_12*cos(x(2)-x(1)); 
       a_12*cos(x(2)-x(1)) a_22];
  F = [0 -a_12*sin(x(2)-x(1));
       a_12*sin(x(2)-x(1)) 0];
  B = [-b_1 0; 0 -b_2];
  % change c depending on analysis
  c = [1;0];
  % c = [-1;1];
  % compute derivative dynamics
  dx = zeros(4,1);
  dd_theta = inv(A)*(c*u-F*[x(3)^2;x(4)^2]-B*[sin(x(1));sin(x(2))]);
  dx(1) = x(3);
  dx(2) = x(4);
  dx(3) = dd_theta(1);
  dx(4) = dd_theta(2);
end
