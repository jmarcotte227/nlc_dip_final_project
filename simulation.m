options = odeset('OutputFcn', @outfun);
[t,y] = ode45(@odefun, [0, 10], [pi; pi; 0; 0], options);
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

% Model simulation
function dxdt=odefun(t,x)
    [A, F, B, C] = matrices(x);
    U = control_outer(x, A, F, B, C);

    rhs = C*U - F*[x(3)^2; x(4)^2] - B*[sin(x(1)); sin(x(2))];
    accel = A\rhs;
    dxdt = [x(3); x(4); accel(1); accel(2)];
end

% Additionally save the input
function status=outfun(t,x,flag)
    persistent u % declare as static variable
    if strcmp(flag, 'done')
        assignin('base', 'u', u); % save to base workspace
        return;
    end
    if strcmp(flag, 'init')
        [A, F, B, C] = matrices(x);
        U = control_outer(x, A, F, B, C);
        u = U';
    else
        for i=1:4 % artifact of using ode45
            [A, F, B, C] = matrices(x(:,i));
            U = control_outer(x(:,i), A, F, B, C);
            u = [u; U'];
        end
    end
    status = 0; % return 0 to continue
end

% Decoupled second order system
function V=control_inner(x)
    k1 = 1; k2 = 1;
    V=-k1*[x(1); x(2)]-k2*[x(3); x(4)];
end

% Wrap the inner controller
function U=control_outer(x, A, F, B, C)
    U=pinv(C)*(B*[sin(x(1)); sin(x(2))]+control_inner(x))
end

% Calculate matrices
function [A, F, B, C] = matrices(x)
    % Model parameters
    m1=1; r1=0.5; l1=1; I1=1/3*m1*l1^2;
    m2=1; r2=0.5; l2=1; I2=1/3*m2*l2^2;
    g = 9.81;
    
    a11=I1+m2*l1^2;
    a12=m2*r2*l1;
    a22=I2;
    b1=(m1*r1+m2*l1)*g;
    b2=m2*r2*g;

    A=[a11, a12*cos(x(2)-x(1));
       a12*cos(x(2)-x(1)), a22];
    
    F=[0, -a12*sin(x(2)-x(1));
       a12*sin(x(2)-x(1)), 0];
    
    B=[-b1, 0; 0, -b2];
    
    C=[-1; 1];
end
