% parameters 
display = true;

% syms a11 a12 a22 x2 x1 x3 x4 b1 b2 c1 c2 
syms x2 x1 x3 x4 
% syms M
syms M1 M2

m_1 = 1;
m_2 = 1;
I_1 = 1;
I_2 = 1;
l = 2;
r_1=1;
r_2=1;

a_11 = I_1+m_2*l^2;
a_12 = m_2*r_2*l;
a_22 = I_2;
b_1 = (m_1*r_1+m_2*l)*g;
b_2 = (m_2*r_2*g);

A = [a11 a12*cos(x2-x1);
     a12*cos(x2-x1) a22];
B = [-b1 0; 0 -b2];
F = [0 -a12*sin(x2-x1);
     a12*sin(x2-x1) 0];
% C = [c1;c2]
C=[1 -1;0 1];
M = [M1;M2];

f = inv(A)*(C*M-F*[x3^2;x4^2]-B*[sin(x1);sin(x2)]);

% f=-inv(A)*(F*[x3^2;x4^2]+B*[sin(x1);sin(x2)]);
% g=-inv(A)*(C*M);

df3dx1 = diff(f(1),x1);
df3dx2 = diff(f(1),x2);
df3dx3 = diff(f(1),x3);
df3dx4 = diff(f(1),x4);

df4dx1 = diff(f(2),x1);
df4dx2 = diff(f(2),x2);
df4dx3 = diff(f(2),x3);
df4dx4 = diff(f(2),x4);

df3dM1 = diff(f(1),M1);
df3dM2 = diff(f(1),M2);

df4dM1 = diff(f(2),M1);
df4dM2 = diff(f(2),M2);

% df3dM = diff(f(1),M);
% df4dM = diff(f(2),M);

if display
  disp("df3/dx1:")
  pretty(simplify(df3dx1))
  disp("df3/dx2:")
  pretty(simplify(df3dx2))
  disp("df3/dx3:")
  pretty(simplify(df3dx3))
  disp("df3/dx4:")
  pretty(simplify(df3dx4))

  disp("df4/dx1:")
  pretty(simplify(df4dx1))
  disp("df4/dx2:")
  pretty(simplify(df4dx2))
  disp("df4/dx3:")
  pretty(simplify(df4dx3))
  disp("df4/dx4:")
  pretty(simplify(df4dx4))
  
  disp("df3/dM1:")
  pretty(simplify(df3dM1))
  disp("df3/dM2:")
  pretty(simplify(df3dM2))
  disp("df4/dM1:")
  pretty(simplify(df4dM1))
  disp("df4/dM2:")
  pretty(simplify(df4dM2))

end


disp("A_tilde: ")
A_tilde = [0 0 1 0;
     0 0 0 1;
     df3dx1 df3dx2 df3dx3 df3dx4;
     df4dx1 df4dx2 df4dx3 df4dx4];
disp("B_tilde: ")
B_tilde = [0 0;
           0 0;
           df3dM1 df3dM2;
           df4dM1 df4dM2];
% B_tilde = [0;
%            0;
%            df3dM;
%            df4dM];


x1=star_x1;
x2=star_x2;
x3=star_x3;
x4=star_x4;
M1=star_M1;
M2=star_M2;

A_tilde_subs = subs(A_tilde);
B_tilde_subs = subs(B_tilde);

disp("Linearized System: ")
pretty(A_tilde_subs)
pretty(B_tilde_subs)
[v,l]=eig(A_tilde_subs);
disp("eigenvalues: ")
pretty(l)

disp("Characteristic Polynomial")
charpoly(A_tilde_subs)

% Routh Hurwitz Table

tab = [charpoly(1) charpoly(3) charpoly(5);
       charpoly(2) charpoly(4) 0;
       0            0          0;
       0            0          0];

