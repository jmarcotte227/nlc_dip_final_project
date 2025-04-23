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
K = [k1, 0;
    0, k2];


%% Storage Function
syms v1 v2
V = 0.5*[x3,x4]*A*[x3;x4]+0.5*[x1,x2]*K*[x1;x2];
dot_A = [ 0, -a12*sin(x2-x1)*(x4-x3);-a12*sin(x2-x1)*(x4-x3), 0];
dot_V = [x3,x4]*[v1;v2]-[x3,x4]*F*[x3^2;x4^2]+0.5*[x3,x4]*dot_A*[x3;x4];


% simplified energy
% V = 0.5*a11*x3^2+0.5*a22*x4^2;
% dot_V = a11*x3*f(1)+a22*x4*f(2);
% 
disp("Two joint control")
pretty(simplify(subs(dot_V)))
disp("------------------")




%% uncertainty analysis
syms m1 m2 l1 l2 r1 r2 g I1 I2 ddt1 ddt2;
a11 = 1/3*m1*l1^2 + m2*l1^2;
a12 = m2*r2*l1;
b1 = (m1*r1+m2*l1)*g;
b2 = m2*r2*g;
a22 =1/3*m2*l2^2; 
ddt = inv(A)*(C*u-F*[x3^2;x4^2]-B*[sin(x1);sin(x2)]);
unc_dynamics = subs((A*[ddt1;ddt2]+F*[x3^2;x4^2]+B*[sin(x1);sin(x2)]));

%unc_dynamics = subs(ddt);
% linear parameterization of the dynamics wrt mass
disp("Linear parameterization")
pretty(simplify(collect(unc_dynamics,[m1,m2])))
%pretty(simplify(collect(subs(unc_dynamics,[l1, l2, r1,r2, g, x1, x2, x3, x4],[1, 1, 0.5, 0.5, 9.81, 1, 1, 1, 1]), [m1,m2])))