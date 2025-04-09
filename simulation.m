T = 10;
[t,y] = ode45(@dip_dynamics, [0,T], [0,0,0,0]);

figure;
plot(t,y);
