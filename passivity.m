clear all
% parameters 
display = true;

syms a11 a12 a22 x2 x1 x3 x4 b1 b2 c1 c2 
syms x2 x1 x3 x4 k1 k2
% syms M
syms u1 u2

A = [a11 a12*cos(x2-x1);
     a12*cos(x2-x1) a22];
B = [-b1 0; 0 -b2];
F = [0 -a12*sin(x2-x1);
     a12*sin(x2-x1) 0];
% C = [c1;c2]
C=[1 -1;0 1];
u = [u1;u2];

% D = 2*eye(2);
D = zeros(2);


y = [x1;x2];

u'*y;

% f = -inv(A)*(F*[x3^2;x4^2]+B*[sin(x1);sin(x2)])
% g = inv(A)*C*u

% define a storage function
f = inv(A)*(C*u-F*[x3^2;x4^2]-B*[sin(x1);sin(x2)]-D*[x3;x4]);

% energy of the system
V = 0.5*a11*x3^2+a12*x3*x4*cos(x1-x2)+0.5*a22*x4^2+b1*cos(x1)+b2*cos(x2);
% dot_V = a11*x3*f(1)+a12*(f(1)*x4*cos(x1-x2)-2*x3*x4*sin(x1-x2)+f(2)*x3*cos(x1-x2))+a22*x4*f(2)-b1*sin(x1)*x3-b2*sin(x2)*x4;
dot_T = a11*x3*f(1)+a12*(f(1)*x4*cos(x1-x2)-x3^2*x4*sin(x1-x2)+x3*x4^2*sin(x1-x2)+f(2)*x3*cos(x1-x2))+a22*x4*f(2)

dot_V = dot_T+x1*x3+x2*x4;


syms v1 v2
dot_A = [ 0, -a12*sin(x2-x1)*(x4-x3);-a12*sin(x2-x1)*(x4-x3), 0];
dot_V = [x3,x4]*[v1;v2]-[x3,x4]*F*[x3^2;x4^2]+0.5*[x3,x4]*dot_A*[x3;x4]


% simplified energy
% V = 0.5*a11*x3^2+0.5*a22*x4^2;
% dot_V = a11*x3*f(1)+a22*x4*f(2);
% 
disp("Before Controller")
pretty(simplify(dot_V))

% Substitute controller into function
% u1=x4*2*a12*sin(x1-x2)-k1*x3
% u2=x3*x4*a12*sin(x1-x2)-k2*(x4-x3)
% disp("After Controller")
% pretty(simplify(subs(dot_V)))

% generic candidate
% syms alph beta gamma delt
% alph=0.5;
% beta=0.5;
% gamma=4;
% delt=4;
% dot_V = x1*x3+x2*x4+x3*f(1)+x4*f(2);
% pretty(simplify(dot_V))

% robot motion candidate
% V = 0.5*[x3, x4]*mod_A*[x3;x4]+0.5*[x1,x2]*[x1;x2];
% dot_A = [ 0, -a12*sin(x2-x1)*(x4-x3);-a12*sin(x2-x1)*(x4-x3), 0];
% dot_V = [x3, x4]*A*f+0.5*[x3,x4]*mod_dot_A*[x3;x4]+[x1,x2]*[x3;x4];


% V = 0.5*[x3, x4]*mod_A*[x3;x4]+0.5*[x1,x2]*[x1;x2];
% mod_A = [a11, 0; 0, a22];
% mod_dot_A = [0,0;0,0]
% dot_A = [ 0, -a12*sin(x2-x1)*(x4-x3);-a12*sin(x2-x1)*(x4-x3), 0];
% dot_V = [x3, x4]*mod_A*f+0.5*[x3,x4]*mod_dot_A*[x3;x4]+[x1,x2]*[x3;x4];
% disp("simplified")
% pretty(simplify(dot_V))
