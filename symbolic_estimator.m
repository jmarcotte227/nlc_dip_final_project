% Symbolic derivation for linearized control matrices

% Model parameters
m1=1; r1=0.5; l1=1; I1=1/3*m1*l1^2;
m2=1; r2=0.5; l2=1; I2=1/3*m2*l2^2;
g = 9.81;

a11=I1+m2*l1^2;
a12=m2*r2*l1;
a22=I2;
b1=(m1*r1+m2*l1)*g;
b2=m2*r2*g;

syms p [4 1]
syms v [2 1]

A = [a11 a12*cos(p(2)-p(1));
     a12*cos(p(2)-p(1)) a22];
B = [-b1 0; 0 -b2];
F = [0 -a12*sin(p(2)-p(1));
     a12*sin(p(2)-p(1)) 0];
C = [1 0; -1 1];

f = [p(3:4); -A\(F*p(3:4)+B*sin(p(1:2)))];
g = [zeros(2); A\C];
h = f+g*v;
dhdx = jacobian(h, p);
dhdu = jacobian(h, v);

A_ = double(subs(dhdx, p, [0; 0; 0; 0]));
B_ = double(subs(dhdu, p, [0; 0; 0; 0]));
C_ = eye(4);
K_ = place(A_, B_, [-0.6 -0.8 -2 -2]);
L_ = place(A_, C_', [-0.6 -0.8 -2 -2]*5);
save('new.mat', 'dhdx', 'dhdu', 'A_', 'B_', 'C_', 'K_', 'L_');
