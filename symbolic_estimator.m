% Symbolic derivation for linearized control matrices

syms a11 a12 a22 b1 b2 x1 x2 x3 x4 u1 u2
x = [x1; x2; x3; x4];
u = [u1; u2];

A = [a11 a12*cos(x(2)-x(1));
     a12*cos(x(2)-x(1)) a22];
B = [-b1 0; 0 -b2];
F = [0 -a12*sin(x(2)-x(1));
     a12*sin(x(2)-x(1)) 0];
C = [1 0; -1 1];

f = [x(3); x(4); -A\(F*[x(3); x(4)]+B*[sin(x(1)); sin(x(2))])];
g = [0; 0; A\(C*u)];

dfdx = jacobian(f, x)+jacobian(g, x);
dfdu = jacobian(f, u)+jacobian(g, u);

dfdx_subs = subs(dfdx, x, [0; 0; 0; 0]);
dfdu_subs = subs(dfdu, x, [0; 0; 0; 0]);

syms k1 k2 k3 k4 k5 k6 k7 k8
K = [k1 k2 k3 k4; k5 k6 k7 k8];
u = C\(-A*K*x+F*[x(3); x(4)]+B*[sin(x(1)); sin(x(2))]);
dudx = jacobian(u, x);
dudx_subs = -subs(dudx, x, [0; 0; 0; 0]);
save('new.mat', 'dfdx_subs', 'dfdu_subs', 'dudx_subs');  
