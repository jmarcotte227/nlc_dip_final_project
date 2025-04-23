load('new.mat')
options = odeset('OutputFcn', @outfun);
[t,y] = ode45(@odefun, [0 20], [1e-1*ones(2,1); zeros(2,1); zeros(4,1); zeros(2,1)], options);
% [t,y,u]=simulation([0, 10], [1e-1*ones(2,1); zeros(2,1); zeros(4,1); zeros(2,1)]);
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
    y = zeros(length(t),length(x));
    u = zeros(length(t),2);
    y(1,:) = x;
    persistent A_ B_ C_ K_ L_
    if isempty(A_)
        A_ = evalin('base', 'A_');
        B_ = evalin('base', 'B_');
        C_ = evalin('base', 'C_');
        K_ = evalin('base', 'K_');
        L_ = evalin('base', 'L_');
    end
    for i=2:length(t)
        x_ = x(1:4); % true state
        e_ = x(5:8); % error state
        h_ = x_-e_; % estimated state
        y_ = x_+[1e-3*ones(2,1); zeros(2,1)]+1e-3*randn(4,1); % measured state
        s_ = x(9:10); % integral state
        [Ax, Fx, Bx, Cx] = matrices(x_);
        [Ah, Fh, Bh, Ch] = matrices(h_);
        U = control_outer(h_, s_, Ah, Fh, Bh, Ch);
        dxdt = [x_(3:4); Ax\(Cx*U - Fx*x_(3:4).^2 - Bx*sin(x_(1:2)))]; % true state
        dhdt = (A_-B_*K_)*h_; % estimated state
        x = x+dt*[dxdt; dxdt-dhdt-L_*(y_-C_*h_); h_(1:2)];
        y(i,:)=x;
        u(i,:)=U;
    end
end

% Model simulation
function dxdt=odefun(t,x)
    persistent A_ B_ C_ K_ L_
    if isempty(A_)
        A_ = evalin('base', 'A_');
        B_ = evalin('base', 'B_');
        C_ = evalin('base', 'C_');
        K_ = evalin('base', 'K_');
        L_ = evalin('base', 'L_');
    end
    x_ = x(1:4); % true state
    e_ = x(5:8); % error state
    h_ = x_-e_; % estimated state
    y_ = x_+[1e-3*ones(2,1); zeros(2,1)]+1e-6*randn(4,1); % measured state
    s_ = x(9:10); % integral state
    [Ax, Fx, Bx, Cx] = matrices(x_);
    [Ah, Fh, Bh, Ch] = matrices(h_);
    U = control_outer(h_, s_, Ah, Fh, Bh, Ch);
    dxdt1 = [x_(3:4); Ax\(Cx*U - Fx*x_(3:4).^2 - Bx*sin(x_(1:2)))]; % true state
    dhdt2 = (A_-B_*K_)*h_; % estimated state
    dxdt = [dxdt1; dxdt1-dhdt2-L_(:,1:4)*(y_-C_(:,1:4)*h_); h_(1:2)];
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
        U = control_outer(x(1:4), x(9:10), A, F, B, C);
        u = U';
    else
        for i=1:4 % artifact of using ode45
            [A, F, B, C] = matrices(x(1:4,i));
            U = control_outer(x(1:4,i), x(9:10,i), A, F, B, C);
            u = [u; U'];
        end
    end
    status = 0; % return 0 to continue
end

% Decoupled second order system
function V=control_inner(x, s)
    persistent K % declare as static variable
    if isempty(K)
        % Pole placement
        A = [zeros(2), eye(2), zeros(2);
             zeros(2), zeros(2), zeros(2);
             eye(2), zeros(2), zeros(2)];
        B = [zeros(2); eye(2); zeros(2)];
        K = place(A, B, [-0.6 -0.8 -1 -1 -0.5 -0.5]);
    end
    
    V = -K*[x; s];
end

% Wrap the inner controller
function U=control_outer(x, s, A, F, B, C)
    U=C\A*(control_inner(x, s) + A\(F*x(3:4).^2+B*sin(x(1:2))));
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
