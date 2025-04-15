m1=1; r1=0.5; l1=1; I1=1/3*m1*l1^2;
m2=1; r2=0.5; l2=1; I2=1/3*m2*l2^2;
g = 9.81;

a11=I1+m2*l1;
a12=m2*r2*l1;
a22=I2;
b1=(m1*r1+m2*l1)*g;
b2=m2*r2*g;

[t,y] = ode45(@(t,x) odefun(t,x,a11,a12,a22,b1,b2), [0, 100], [0; 0; 0; 0]);
plot(t,y(:,1:2))

function dxdt=odefun(t,x,a11,a12,a22,b1,b2)

A=[a11, a12*cos(x(2)-x(1));
   a12*cos(x(2)-x(1)), a22];

F=[0, -a12*sin(x(2)-x(1));
   a12*sin(x(2)-x(1)), 0];

B=[-b1, 0; 0, -b2];

C=[1, 0; -1, 1];

M=[0; 0];


% rhs = C*M - F*[x(3)^2; x(4)^2] - B*[sin(x(1)); sin(x(2))];
% accel = inv(A)*rhs;
% with damping
D = 0.1*eye(2); % must be positive semi-definite
rhs = C*M - F*[x(3)^2; x(4)^2] - B*[sin(x(1)); sin(x(2))]-D*[x(3);x(4)];
accel = inv(A)*rhs;
dxdt = [x(3); x(4); accel(1); accel(2)];
end
