load('new.mat')
% options = odeset('OutputFcn', @outfun,'RelTol', 1e-3, 'AbsTol', 1e-4);
% [t,y] = ode23(@odefun, 0:0.01:5, [pi; pi; 0; 0], options);
[t,y,u]=simulation([0, 10], [1e-1; 1e-1; 0; 0]);
figure(1)
plot(t,y(:,1:2))
xlabel('t (s)')
ylabel('\theta (rad)')
title('State Variables')
legend('\theta_1', '\theta_2')
figure(2)
plot(t,u)
xlabel('t (s)')
ylabel('u (N*m)')
title('Control Effort')
legend('u_1', 'u_2')

% Fixed time step simulation
function [t,y,u]=simulation(t,x)
    dt = 0.001;
    t = t(1):dt:t(end);
    y = zeros(length(t),4);
    u = zeros(length(t),2);
    y(1,:) = x;
    for i=2:length(t)
        [A, F, B, C] = matrices(x);
        U = control_outer(x, A, F, B, C);
        dxdt = [x(3); x(4); A\(C*U - F*x(3:4).^2 - B*sin(x(1:2)))];
        x = x+dt*dxdt;
        y(i,:)=x;
        u(i,:)=U;
    end
end

% Model simulation
function dxdt=odefun(t,x)
    [A, F, B, C] = matrices(x);
    U = control_outer(x, A, F, B, C);
    dxdt = [x(3); x(4); A\(C*U - F*x(3:4).^2 - B*sin(x(1:2)))];
    t
end

% Additionally save the input
function status=outfun(t,x,flag)
    persistent u % declare as static variable
    if strcmp(flag, 'done')
        assignin('base', 'u', u); % save to base workspace
        return;
    end
    if strcmp(flag, 'init')
        [A, F, B, C] = matrices(x(1:4));
        U = control_outer(x(1:4), A, F, B, C);
        u = U';
    else
        [A, F, B, C] = matrices(x);
        U = control_outer(x, A, F, B, C);
        u = [u; U'];
    end
    status = 0; % return 0 to continue
end

% Decoupled second order system
function V=control_inner(x)
    persistent A B K L C
    if isempty(A)
        A = evalin('base', 'dfdx_subs');
        B = evalin('base', 'dfdu_subs');
        K = evalin('base', 'dudx_subs');
    
        % Model parameters
        m1=1; r1=0.5; l1=1; I1=1/3*m1*l1^2;
        m2=1; r2=0.5; l2=1; I2=1/3*m2*l2^2;
        g = 9.81;
        
        a11_=I1+m2*l1^2;
        a12_=m2*r2*l1;
        a22_=I2;
        b1_=(m1*r1+m2*l1)*g;
        b2_=m2*r2*g;
    
        % LQR control
        A = [0 0 1 0;
             0 0 0 1;
             0 0 0 0;
             0 0 0 0];
         
        B = [0 0;
             0 0;
             1 0;
             0 1];
        
        Q = diag([5, 1, 1, 1]);
        R = diag([10, 1]);
        K_ = lqr(A, B, Q, R);
        k1_ = K_(1,1); k2_ = K_(1,2); k3_ = K_(1,3); k4_ = K_(1,4);
        k5_ = K_(2,1); k6_ = K_(2,2); k7_ = K_(2,3); k8_ = K_(2,4);
    
        syms a11 a12 a22 b1 b2 k1 k2 k3 k4 k5 k6 k7 k8
        A = double(subs(A, [a11; a12; a22; b1; b2], [a11_; a12_; a22_; b1_; b2_]));
        B = double(subs(B, [a11; a12; a22], [a11_; a12_; a22_]));
        K = double(subs(K, [a11; a12; a22; b1; b2; k1; k2; k3; k4; k5; k6; k7; k8], [a11_; a12_; a22_; b1_; b2_; k1_; k2_; k3_; k4_; k5_; k6_; k7_; k8_]));
        L = A+diag([2,1,1,1]);
        C = eye(4);
    end

    V = -K*x;
end

% Wrap the inner controller
function U=control_outer(x, A, F, B, C)
    x = x+1e-3*randn(4,1);
    U=C\A*(control_inner(x) + A\(F*x(3:4).^2+B*sin(x(1:2))));
end

% Calculate matrices
function [A, F, B, C] = matrices(x)
    persistent m1 r1 l1 I1 m2 r2 l2 I2 g a11 a12 a22 b1 b2
    if isempty(g)
        % Model parameters
        m1=1; r1=0.5; l1=1; I1=1/3*m1*l1^2;
        m2=1; r2=0.5; l2=1; I2=1/3*m2*l2^2;
        g = 9.81;
        
        a11=I1+m2*l1^2;
        a12=m2*r2*l1;
        a22=I2;
        b1=(m1*r1+m2*l1)*g;
        b2=m2*r2*g;
    end

    A=[a11, a12*cos(x(2)-x(1));
       a12*cos(x(2)-x(1)), a22];
    
    F=[0, -a12*sin(x(2)-x(1));
       a12*sin(x(2)-x(1)), 0];
    
    B=[-b1, 0; 0, -b2];
    
    C=[1, -1; 0, 1];
end
