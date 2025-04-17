global Df_1
global Df_2
global hm1
global hm2
global sdxdt
global sx
global sU

syms hm1 hm2 
sdxdt = sym("dx", [6 1]);
sx = sym("x", [8 1]);
sU = sym("U", [2 1]);
est_e = f_m(sx, sdxdt, sU, hm1, hm2);
Df_1 = diff(est_e'*est_e,hm1);
Df_2 = diff(est_e'*est_e,hm2);




function [t,y,u]=simulation(t,x)
    dt = 0.01;
    t = t(1):dt:t(end);
    y = zeros(length(t),8);
    u = zeros(length(t),2);
    y(1,:) = x;
    for i=2:length(t)
        global Df_0
        global Df_1
        global hm0
        global hm1
        global sdxdt
        global sx
        global sU
        % actual dynamics
        m1 = 1;
        m2 = 1;
        p=0.001;
        [A, F, B, C] = matrices(x, m1, m2);
        % uncertain dynamics
        % uncertain mass
        hat_m1 = x(5);
        hat_m2 = x(6);
        hB = unc_B(hat_m1, hat_m2);
    
        U = passive_control(x, C, hB);
    
        % Integral Control
        % U = int_passive_control(x, C, hat_B);
    
        % U = control_outer(x, A, F, B, C);
    
        rhs = C*U - F*[x(3)^2; x(4)^2] - B*[sin(x(1)); sin(x(2))];
        accel = A\rhs;
        dxdt = [x(3); x(4); accel(1); accel(2); 0; 0];
        
        % substitute 
        grad_m1 = subs(Df_1, [hm1, hm2], [hat_m1, hat_m2]);
        grad_m1 = subs(grad_m1, [sdxdt, sx], [dxdt, x]);
        grad_m1 = subs(grad_m1, sU, U);
        grad_m2 = subs(Df_2, [hm1, hm2], [hat_m1, hat_m2]);
        grad_m2 = subs(grad_m2, [sdxdt, sx], [dxdt, x]);
        grad_m2 = subs(grad_m2, sU, U);

    % update parameters
        dxdt(7) = p*grad_m1;
        dxdt(8) = p*grad_m2;
        % dxdt = [x(3); x(4); accel(1); accel(2); 0; 0
        x = x+dt*dxdt;
        y(i,:) = x;
        u(i,:) = U;
        disp(t(i))
    end
end 
% options = odeset('OutputFcn', @outfun);
% [t,y] = ode45(@odefun, [0, 5], [-pi; -pi; 0; 0; 1; 1]);
[t, y, u] = simulation([0, 5], [-pi; -pi; 0; 0; 1.1; 1.1]);
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
figure(3)
% nominal B
B = unc_B(1,1);
plot(t,y(:,5:6))
xlabel('t (s)')
ylabel('m (kg)')
title('Parameter update')
legend('b_1', 'b_2')

% Model simulation
% function dxdt=odefun(t,x)
%     % actual dynamics
%     m1 = 1;
%     m2 = 1;
%     p=0.5;
%     global Df_1
%     global Df_2
%     global hm1
%     global hm2
%     global sdxdt
%     global sx
%     global sU
%     [A, F, B, C] = matrices(x, m1, m2);
%     % uncertain dynamics
%     % uncertain mass
%     hat_m1 = 1;
%     hat_m2 = 1;
%     B
%     hB = unc_B(hat_m1, hat_m2)
% 
%     U = passive_control(x, C, hB)
% 
%     % Integral Control
%     % U = int_passive_control(x, C, hat_B);
% 
%     % U = control_outer(x, A, F, B, C);
% 
%     rhs = C*U - F*[x(3)^2; x(4)^2] - B*[sin(x(1)); sin(x(2))];
%     accel = A\rhs;
%     dxdt = [x(3); x(4); accel(1); accel(2); 1; 1];
% 
%     % substitute 
% %     grad_m1 = subs(Df_1, [hm1, hm2], [hat_m1, hat_m2]);
% %     grad_m1 = subs(grad_m1, [sdxdt, sx], [dxdt, x]);
% %     grad_m1 = subs(grad_m1, sU, U);
% %     grad_m2 = subs(Df_2, [hm1, hm2], [hat_m1, hat_m2]);
% %     grad_m2 = subs(grad_m2, [sdxdt, sx], [dxdt, x]);
% %     grad_m2 = subs(grad_m2, sU, U);
% % 
% % % update parameters
% %     dxdt(5) = p*grad_m1;
% %     dxdt(6) = p*grad_m2;
%     % dxdt = [x(3); x(4); accel(1); accel(2); 0; 0
% end

% compute the error between the estimated and real system
function est_e = f_m(dt,x,dxdt, U, hm1, hm2)
    % compute estimated dynamics
    [hA, hF, hB, hC] = matrices(x, hm1, hm2);
    rhs = hC*U - hF*[x(3)^2; x(4)^2] - hB*[sin(x(1)); sin(x(2))];
    accel = hA\rhs;
    dhxdt = [x(3); x(4); accel(1); accel(2); -x(1); -x(2)];
    est_e = (x+(dt*dxdt)) - (x_hat)dhxdt;
end
% Additionally save the input
% function status=outfun(t,x,flag)
%     persistent u % declare as static variable
%     if strcmp(flag, 'done')
%         assignin('base', 'u', u); % save to base workspace
%         return;
%     end
%     if strcmp(flag, 'init')
%         [A, F, B, C] = matrices(x);
%         U = control_outer(x, A, F, B, C);
%         u = U';
%     else
%         for i=1:4 % artifact of using ode45
%             [A, F, B, C] = matrices(x(:,i));
%             U = control_outer(x(:,i), A, F, B, C);
%             u = [u; U'];
%         end
%     end
%     status = 0; % return 0 to continue
% end

% Decoupled second order system
function V=control_inner(x)
    zeta1=1;
    omega1=1;
    zeta2=1;
    omega2=1;
    
    k11=omega1^2;
    k21=2*zeta1*omega1;
    k12=omega2^2;
    k22=2*zeta2*omega2;
    
    V = [-k11*x(1)-k21*x(3); -k12*x(2)-k22*x(4)];
end

% Wrap the inner controller
function U=control_outer(x, A, F, B, C)
    U=C\A*(control_inner(x) + A\(F*[x(3)^2; x(4)^2]+B*[sin(x(1)); sin(x(2))]));
end

% Passivity-Based control
function U=passive_control(x, C, B)
    % gains
    K1 = 2*eye(2);
    K2 = 2*eye(2);
    U=C\(B*[sin(x(1)); sin(x(2))]-K1*[x(1);x(2)]-K2*[x(3);x(4)]);
end
% Passivity-Based control with integrator
function U=int_passive_control(x, C, B)
    % gains
    K1 = eye(2);
    K2 = eye(2);
    K3 = 0.5*eye(2);
    U=C\(B*[sin(x(1)); sin(x(2))]-K1*[x(1);x(2)]-K2*[x(3);x(4)]+K3*[x(5);x(6)]);
end

% Calculate matrices
function [A, F, B, C] = matrices(x, m1, m2)
    % Model parameters
    r1=0.5; l1=1; I1=1/3*m1*l1^2;
    r2=0.5; l2=1; I2=1/3*m2*l2^2;
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
    
    C=[1, -1; 0, 1];
end

% calculate uncertain b
function hat_B = unc_B(hat_m1, hat_m2)
    % known parameters
    r1=0.5; l1=1; r2=0.5;
    g=9.81;

    b1=(hat_m1*r1+hat_m2*l1)*g;
    b2=hat_m2*r2*g;
    hat_B=[-b1, 0; 0, -b2];
end
