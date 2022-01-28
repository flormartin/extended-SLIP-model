clc; clear all; close all;
addpath('C:\Program Files\MATLAB\R2021a\toolbox\casadi-windows-matlabR2016a-v355')
import casadi.*

%% Initialization

opti = casadi.Opti();

% decision_variables = [X, T_flight, T_stance]

% Number of nodes
N = 80;
N_phase = N/2;

% x = [x y phi theta dx dy dphi dtheta]

X = opti.variable(8, N+1);
T_st = opti.variable(1);
T_fl = opti.variable(1);
% actuation parameters
phi0 = opti.variable(1);
alpha0 = opti.variable(1);
omega = opti.variable(1);

% Parameters
global L0 l0 l1 M m_B m_L g c_phi d1 c_theta d2 c1 t_apex;
L0=1; M=80; g=9.81;
l0 = .6*L0;                     %lower leg length
l1 = .4*L0;                     %upper leg length
l_l = l0;   % lower leg
l_u = l1;   % upper leg
m_L = .32*M;                    %leg mass
m_B = M-m_L;                    %rest of body mass
c_phi = 750;                    %Nm/rad
d1 = 2*sqrt(c_phi*m_L);
c_theta = 1000;
d2 = 2*sqrt(c_theta*m_L);
c1 = 20 * 1000;

c1 = 20*1000;
%alpha0 = (90 - 62) * pi/ 180.;
t_apex = 0;
%omega = 50*pi/180; %rad/s

%% Flight phase - Mode 2

h_fl = T_fl/N_phase;    % step width

% Runge-Kutta
for k = 1:N_phase
    
    t = k * h_fl;
    k1 = mode2Florian(t, X(:,k), alpha0, omega);
    k2 = mode2Florian(t, X(:,k) + 0.5 * h_fl * k1, alpha0, omega);
    k3 = mode2Florian(t, X(:,k) + 0.5 * h_fl * k2, alpha0, omega);
    k4 = mode2Florian(t, X(:,k) + h_fl * k3, alpha0, omega);
    
    x_next = X(:,k) + 1./6. * h_fl * (k1 + 2 * k2 + 2 * k3 + k4);
    
    opti.subject_to(X(:,k+1) == x_next);
end

phi = X(3,N_phase+1);
theta = X(4,N_phase+1);
offset = X(1,N_phase+1) + l1 * sin(phi) + l0 * sin(phi - theta);

%% Stance phase - Mode 1

h_st = T_st/N_phase;    % step width

% Runge-Kutta
for k = N_phase+1:N

    k1 = mode1Florian(X(:,k), phi0, offset);
    k2 = mode1Florian(X(:,k) + 0.5 * h_st * k1, phi0, offset);
    k3 = mode1Florian(X(:,k) + 0.5 * h_st * k2, phi0, offset);
    k4 = mode1Florian(X(:,k) + h_st * k3, phi0, offset);
    
    x_next = X(:,k) + 1./6. * h_st * (k1 + 2 * k2 + 2 * k3 + k4);
    
    opti.subject_to(X(:,k+1) == x_next);
end

%% Periodicity conditions

% x = [x y phi theta dx dy dphi dtheta]
opti.subject_to(X(2,1) == X(2,N+1));
opti.subject_to(X(3,1) == X(3,N+1));
opti.subject_to(X(4,1) == X(4,N+1));
opti.subject_to(X(5,1) == X(5,N+1));
opti.subject_to(X(6,1) == X(6,N+1));
opti.subject_to(X(7,1) == X(7,N+1));
opti.subject_to(X(8,1) == X(8,N+1));


%% Guard conditions
% l0: length lower leg; l1: length upper leg

% liftoff
x_lo_knee = (X(1,N+1) - offset) + l1 * sin(X(3,N+1)); 
y_lo_knee = X(2,N+1) - l1 * cos(X(3,N+1)) ;
% constraint: length of lower leg equal to 0.6
opti.subject_to(sqrt(x_lo_knee^2 + y_lo_knee^2) == l0);  

% touchdown  
y_td_knee = X(2,N_phase+1) - l1 * cos(X(3,N_phase+1));   
alpha_td_knee = X(3,N_phase+1) - X(4,N_phase+1);%%%%%%%%%%%%%%%%%%
% constraint: height of knee computed from hip = height of knee computed from lower leg 
opti.subject_to(y_td_knee == l0 * cos(alpha_td_knee));

% constraint: height of foot >= 0 at flight phase
for i = 1:N_phase
    opti.subject_to( (X(2,i) - l1 * cos(X(3,i)) - cos(X(3,i)-X(4,i)) * l0) >= 0);
end

% l_s - length of spring (lower leg)
% actually l_s should equal to 0.6 since at N_phase+1 it hasn't shrunk yet
l_s = sqrt((offset-(X(1,N_phase+1)+l1*sin(X(3,N_phase+1))))^2+(X(2,N_phase+1)-l1*cos(X(3,N_phase+1)))^2);

% foot on the ground at N_phase+1 (y_foot = 0)
opti.subject_to( (X(2,N_phase+1) - l1 * cos(X(3,N_phase+1)) - cos(X(3,N_phase+1)-X(4,N_phase+1)) * l_s) == 0);
% foot start from the ground at time 1 (y_foot(time 1) = 0)
opti.subject_to( (X(2,1) - l1 * cos(X(3,1)) - cos(X(3,1)-X(4,1)) * l_s) == 0);


% %% Path contraints
% % x = [x y phi theta dx dy dphi dtheta]
% opti.subject_to(T_fl >= 0.15);  
% opti.subject_to(T_st >= 0.15);
% 
% opti.subject_to(X(1,1) == 0); % starting point - horizontal coordinate x
% opti.subject_to(X(2,:) >= .4); % height of hip
% opti.subject_to(X(3,:) >= -pi*3/4);   % phi 
% opti.subject_to(X(3,:) <= pi*3/4);
% opti.subject_to(X(3,1) <= 0);
% opti.subject_to(X(4,:) >= 0);
% opti.subject_to(X(4,1) == 0);
% opti.subject_to(X(4,N_phase) == 0);
% opti.subject_to(X(5,:) >= 0); % dx - hip keeps moving forward not backward
% 
% % new constraints including the actuation parameters as optimization
% % variables
% opti.subject_to(phi0 <= -pi/16);
% opti.subject_to(phi0 >= -pi/2);
% opti.subject_to(alpha0 <= pi/2);
% opti.subject_to(alpha0 >= -pi/2);
% opti.subject_to(omega >= 10*pi/180);

%% Initial guess

% for relaxed solution - without guard 
x0 = [0.; .8; -35 * pi/ 180.; 0.; 5; .5; 0; 0];
opti.set_initial(X, repmat(x0,1,N+1));
opti.set_initial(phi0, -35*pi/180);
opti.set_initial(alpha0, -35*pi/180);
opti.set_initial(omega, 50*pi/180);
opti.set_initial(T_fl, 0.2);
opti.set_initial(T_st, 0.2);

% for solution with guard

% load('relaxed_solution_2.mat')
% load('strict_solution_2.mat')
% opti.set_initial(X, x_sol);
% opti.set_initial(T_fl, time_fl);
% opti.set_initial(T_st, time_st);

% show progress of optimization 
opti.callback(@(i) plot(opti.debug.value(X(1,1:N+1)), opti.debug.value(X(2,1:N+1)),'b') )

opti.solver('ipopt');   % interior point method
sol = opti.solve();

%% Test & relaxed solution
% 
% Y_td_knee = sol.value(X(2,N_phase+1)) - l_u * cos(sol.value(X(3,N_phase+1))) % touch down    
% Angle_td_knee = sol.value(X(3,N_phase+1)) - sol.value(X(4,N_phase+1));%%%%%%%%%
% l_l * cos(Angle_td_knee)
% 
% % Offset = sol.value(X(1,N_phase+1)) + L0*sin(sol.value(X(3,N_phase+1)));
% Offset = sol.value(X(1,N_phase+1)) + l_u * sin(sol.value(X(3,N_phase+1))) + l_l * sin(sol.value(Angle_td_knee));
% X_lo_knee = (sol.value(X(1,N+1)) - Offset) + l_u * sin(sol.value(X(3,N+1)));
% Y_lo_knee = sol.value(X(2,N+1)) - l_u * cos(sol.value(X(3,N+1))) ;
% sqrt(X_lo_knee^2 + Y_lo_knee^2)
% % 
% % sol.value(offset)
% 
% time_fl = sol.value(T_fl)
% time_st = sol.value(T_st)
% x_sol = sol.value(X)

%% Plot 

x_hip_sol = sol.value(X(1,:));
y_hip_sol = sol.value(X(2,:));
phi_sol = sol.value(X(3,:));
theta_sol = sol.value(X(4,:));
angle_sol = phi_sol - theta_sol;%%%%%%%%%%
phi_td = sol.value(phi_sol(N_phase+1)); 
theta_td = sol.value(theta_sol(N_phase+1)); 
angle_td = sol.value(phi - theta);%%%%%%%%%%%
offset_sol = x_hip_sol(N_phase+1) + l_u * sin(phi_td) + l_l * sin(angle_td);

% Knee
x_knee = zeros(1,N+1);
y_knee = zeros(1,N+1);
x_knee = x_hip_sol + l_u * sin(phi_sol);
y_knee = y_hip_sol - l_u * cos(phi_sol);

% Foot 
x_foot = zeros(1,N+1);
y_foot = zeros(1,N+1);
l_s_sol = zeros(1,N+1);
l_s_sol(1:N_phase) = 0.6;
for i = N_phase+1:N+1
    l_s_sol(i) = sqrt((offset_sol-(x_hip_sol(i)+l1*sin(phi_sol(i))))^2+(y_hip_sol(i)-l1*cos(phi_sol(i)))^2);
end
x_foot = x_knee + l_s_sol .* sin(angle_sol);
y_foot = y_knee - l_s_sol .* cos(angle_sol);

 
figure; hold on;
for i = 1:N+1
    plot(x_hip_sol(i), y_hip_sol(i),'o','Color','k');   % position of hip
    plot(x_knee(i), y_knee(i),'o','Color','b'); % position of knee
    plot(x_foot(i), y_foot(i),'o','Color','k');  % position of foot
    plot([x_hip_sol(i),x_knee(i)],[y_hip_sol(i),y_knee(i)],'Color','k');    % upper leg
    plot([x_knee(i),x_foot(i)],[y_knee(i),y_foot(i)],'Color','k');  % lower leg
end
plot(x_foot, y_foot,'-','Color','r');  %  foot
plot(x_knee, y_knee,'-','Color','r');  %  knee
plot(x_hip_sol, y_hip_sol,'-','Color','r');  %  hip

% 
% 
%% Test 
% 
% 
% load('strict_solution_3.mat')
% x_hip_sol = x_sol(1,:);
% y_hip_sol = x_sol(2,:);
% phi_sol = x_sol(3,:);
% theta_sol = x_sol(4,:);
% angle_sol = phi_sol + theta_sol;
% phi_td = phi_sol(N_phase+1); 
% theta_td = theta_sol(N_phase+1); 
% angle_td = phi - theta;%%%%%%%%%%%
% offset_sol = x_hip_sol(N_phase+1) + l_u * sin(phi_td) + l_l * sin(angle_td);
% 
% % Knee
% x_knee = zeros(1,N);
% y_knee = zeros(1,N);
% x_knee = x_hip_sol + l_u * sin(phi_sol);
% y_knee = y_hip_sol - l_u * cos(phi_sol);
% 
% % Foot -
% x_foot = zeros(1,N);
% y_foot = zeros(1,N);
% x_foot(1:N) = x_knee(1:N) + l_l * sin(angle_sol(1:N));
% y_foot(1:N) = y_knee(1:N) - l_l * cos(angle_sol(1:N));
% 
% 
% % 
% figure; hold on;
% for i = 1:N
%     plot(x_hip_sol(i), y_hip_sol(i),'o','Color','k');
%     plot(x_knee(i), y_knee(i),'o','Color','b');
%     plot(x_foot(i), y_foot(i),'o','Color','k');
%     plot([x_hip_sol(i),x_knee(i)],[y_hip_sol(i),y_knee(i)],'Color','k');
%     plot([x_knee(i),x_foot(i)],[y_knee(i),y_foot(i)],'Color','k');
% end
