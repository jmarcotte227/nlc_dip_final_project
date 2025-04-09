T = 100;
[t,y] = ode45(@dip_dynamics, [0,T], [pi,0,0,0]);

figure;
plot(t,y(:,1:2));
