% parameters 
display = true;

syms a11 a12 a22 x2 x1 x3 x4 b1 b2 c1 c2 
syms x2 x1 x3 x4 
% syms M
syms u1 u2

m1=1; r1=0.5; l1=1; I1=1/3*m1*l1^2;
m2=1; r2=0.5; l2=1; I2=1/3*m2*l2^2;
g = 9.81;

a11=I1+m2*l1^2;
a12=m2*r2*l1;
a22=I2;
b1=(m1*r1+m2*l1)*g;
b2=m2*r2*g;

A = [a11 a12*cos(x2-x1);
     a12*cos(x2-x1) a22];
B = [-b1 0; 0 -b2];
F = [0 -a12*sin(x2-x1);
     a12*sin(x2-x1) 0];
% C = [c1;c2]
C=[1 -1;0 1];
M = [u1;u2];

f = inv(A)*(C*M-F*[x3^2;x4^2]-B*[sin(x1);sin(x2)]);
