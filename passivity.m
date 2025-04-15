clear all
% parameters 
display = true;

syms a11 a12 a22 x2 x1 x3 x4 b1 b2 c1 c2 
syms x2 x1 x3 x4 
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


y = [x1;x2];

u'*y;

% f = -inv(A)*(F*[x3^2;x4^2]+B*[sin(x1);sin(x2)])
% g = inv(A)*C*u

% define a storage function
f = inv(A)*(C*u-F*[x3^2;x4^2]-B*[sin(x1);sin(x2)]);

% energy of the system
% V = 0.5*a11*x3^2+a12*x3*x4*cos(x1-x2)+0.5*a22*x4^2+b1*cos(x1)+b2*cos(x2);
% dot_V = a11*x3*f(1)+a12*(f(1)*x4*cos(x1-x2)-2*x3*x4*sin(x1-x2)+f(2)*x3*cos(x1-x2))+a22*x4*f(2)-b1*sin(x1)*x3-b2*sin(x2)*x4;

% simplified energy
V = 0.5*a11*x3^2+0.5*a22*x4^2;
dot_V = a11*x3*f(1)+a22*x4*f(2);
% 
disp("Before Controller")
pretty(collect(simplify(dot_V), [u1,u2]))

% Substitute controller into function
% u1=x4*2*a12*sin(x1-x2)
% u2=x3*x4*a12*sin(x1-x2)
% disp("After Controller")
% pretty(simplify(subs(dot_V)))

% generic candidate
% syms alph beta gamma delt
% alph=0.5;
% beta=0.5;
% gamma=4;
% delt=4;
dot_V = x1*x3+x2*x4+x3*f(1)+x4*f(2);
pretty(simplify(dot_V))
