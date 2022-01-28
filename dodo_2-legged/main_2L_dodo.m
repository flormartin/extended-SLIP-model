clc; clear all; close all;
addpath('C:\Users\loyst\OneDrive - Technische Universitat Munchen\Studium\Master_EPM-in-Maschinenbau_TUM\05_WS21-22\02_Fächer\09_Dodo-Alive-MechanismDesign\Casadi')
import casadi.*


%% Initialization

opti = casadi.Opti();

% decision_variables = [X, T_flight, T_stance]

% Number of nodes
N = 80;
N_phase = N/2;

% x = [x y phi phi1 theta theta1 dx dy dphi dphi1 dtheta dtheta1]

X = opti.variable(12, N+1);
T_st_L1 = opti.variable(1);
T_st_L2 = opti.variable(1);
% actuation parameters
phi0 = opti.variable(1);
theta0 = opti.variable(1);
alpha0 = opti.variable(1);
omega = opti.variable(1);
% control input
U = opti.variable(1, N_phase);  % ddphi

% Parameters
global L0 l0 l1 M m_B m_L g c_hip dhip c_knee dknee c1 offset;
L0=.40; M=12; g=9.81;
l0 = .6966*L0;                     %lower leg length
l1 = .3034*L0;                     %upper leg length
m_L = .18*M;                    %leg mass
m_B = M-m_L;                  %rest of body mass
c_hip = 1125;                    %Nm/rad
dhip = 2*sqrt(c_hip*m_L);
c_knee = 150;
dknee = 2*sqrt(c_knee*m_L);
c1 = 3000;

%alpha0 = (90 - 62) * pi/ 180.;
%t_apex = X(6,1);
%omega = 50*pi/180; %rad/s

%% Stance phase - Mode 1

h_st = T_st_L1/N_phase;    % step width
offset = X(1,0) + l1 * sin(X(3,0)) + l0 * sin(X(3,0) - X(4,0));
% Runge-Kutta
for k = N_phase+1:N
    j = k-N_phase;

    k1 = mode1_2L_dodo(X(:,k), phi0, theta0, offset, U(j));
    k2 = mode1_2L_dodo(X(:,k) + 0.5 * h_st * k1, phi0, theta0, offset, U(j));
    k3 = mode1_2L_dodo(X(:,k) + 0.5 * h_st * k2, phi0, theta0, offset, U(j));
    k4 = mode1_2L_dodo(X(:,k) + h_st * k3, phi0, theta0, offset, U(j));
    
    x_next = X(:,k) + 1./6. * h_st * (k1 + 2 * k2 + 2 * k3 + k4);
    
    opti.subject_to(X(:,k+1) == x_next);
end

theta = X(5,N_phase+1);
theta1 = X(6,N_phase+1);
offset = X(1,N_phase+1) + l1 * sin(theta) + l0 * sin(theta - theta1);

%% Stance phase of leg 2 - Mode 2

h_st = T_st_L2/N_phase;    % step width

% Runge-Kutta
for k = N_phase+1:N
    j = k-N_phase;

    k1 = mode2_2L_dodo(X(:,k),phi0, theta0, offset, U(j));
    k2 = mode2_2L_dodo(X(:,k) + 0.5 * h_st * k1, phi0, theta0, offset, U(j));
    k3 = mode2_2L_dodo(X(:,k) + 0.5 * h_st * k2, phi0,theta0, offset, U(j));
    k4 = mode2_2L_dodo(X(:,k) + h_st * k3, phi0, theta0, offset, U(j));
    
    x_next = X(:,k) + 1./6. * h_st * (k1 + 2 * k2 + 2 * k3 + k4);
    
    opti.subject_to(X(:,k+1) == x_next);
end

%% Periodicity conditions

% x = [x y phi phi1 theta theta1 dx dy dphi dphi1 dtheta dtheta1]
% opti.subject_to(X(2,1) == X(2,N+1));
% opti.subject_to(X(3,1) == X(3,N+1));
% opti.subject_to(X(4,1) == X(4,N+1));
% opti.subject_to(X(5,1) == X(5,N+1));
% opti.subject_to(X(6,1) == X(6,N+1));
% opti.subject_to(X(7,1) == X(7,N+1));
% opti.subject_to(X(8,1) == X(8,N+1));
% opti.subject_to(X(9,1) == X(9,N+1));
% opti.subject_to(X(10,1) == X(10,N+1));
% opti.subject_to(X(11,1) == X(11,N+1));
% opti.subject_to(X(12,1) == X(12,N+1));

%% Guard conditions
% l0: length lower leg; l1: length upper leg

% touchdown leg 2
y_td_knee = X(2,N_phase+1) - l1 * cos(X(5,N_phase+1));   
alpha_td_leg2 = X(5,N_phase+1) - X(6,N_phase+1);%%%%%%%%%%%%%%%%%%
% constraint: height of knee computed from hip = height of knee computed from lower leg 
opti.subject_to(y_td_knee == l0 * cos(alpha_td_leg2));

% touchdown leg 1 
y_td_knee = X(2,N_phase+1) - l1 * cos(X(3,N_phase+1));   
alpha_td_leg1 = X(3,N_phase+1) - X(4,N_phase+1);%%%%%%%%%%%%%%%%%%
% constraint: height of knee computed from hip = height of knee computed from lower leg 
opti.subject_to(y_td_knee == l0 * cos(alpha_td_leg1));

% % constraint: height of foot >= 0 at flight phase
% for i = 1:N_phase
%     opti.subject_to( (X(2,i) - l1 * cos(X(3,i)) - cos(X(3,i)-X(4,i)) * l0) >= 0);
% end

% % l_s - length of spring (lower leg)
% % actually l_s should equal to 0.6 since at N_phase+1 it hasn't shrunk yet
% l_s = sqrt((offset-(X(1,N_phase+1)+l1*sin(X(3,N_phase+1))))^2+(X(2,N_phase+1)-l1*cos(X(3,N_phase+1)))^2);
% 
% % foot on the ground at N_phase+1 (y_foot = 0)
% opti.subject_to( (X(2,N_phase+1) - l1 * cos(X(3,N_phase+1)) - cos(X(3,N_phase+1)-X(4,N_phase+1)) * l_s) == 0);
% % foot start from the ground at time 1 (y_foot(time 1) = 0)
% opti.subject_to( (X(2,1) - l1 * cos(X(3,1)) - cos(X(3,1)-X(4,1)) * l_s) == 0);

% 
%% Path contraints
% x = [x y phi phi1 theta theta1 dx dy dphi dphi1 dtheta dtheta1]
% opti.subject_to(T_st_L1 >= 0.15);  
% opti.subject_to(T_st_L2 >= 0.15);
% 
% opti.subject_to(X(1,1) == 0); % starting point - horizontal coordinate x
% opti.subject_to(X(2,:) >= .4); % height of hip
% opti.subject_to(X(3,:) >= -pi*3/4);   % phi 
% opti.subject_to(X(3,:) <= pi*3/4);
% opti.subject_to(X(3,1) <= 0);
% opti.subject_to(X(4,:) >= 0);
% opti.subject_to(X(4,1) == 0);
% opti.subject_to(X(4,N_phase) == 0);
% opti.subject_to(X(7,:) >= 0); % dx - hip keeps moving forward not backward
% %opti.subject_to(X(8,1) >= 0); % dy, unsure about that one
% 
% % new constraints including the actuation parameters as optimization
% % variables
% opti.subject_to(phi0 <= -pi/16);
% opti.subject_to(phi0 >= -pi/2);
% opti.subject_to(theta0 <= -pi/16);
% opti.subject_to(theta0 >= -pi/2);
% opti.subject_to(alpha0 <= pi/2);
% opti.subject_to(alpha0 >= -pi/2);
% opti.subject_to(omega >= 10*pi/180);

%% Initial guess
% % 
% % for relaxed solution - without guard 
x0 = [0.; .4; 20 * pi/ 180.; 0.; -20 * pi/ 180.; 0.; 2; -.2; 0; 0; 0; 0];
opti.set_initial(X, repmat(x0,1,N+1));
opti.set_initial(phi0, 20*pi/180);
opti.set_initial(theta0, -20*pi/180);
opti.set_initial(alpha0, 20*pi/180);
opti.set_initial(omega, 10*pi/180);
opti.set_initial(T_st_L1, 2);
opti.set_initial(T_st_L2, 2);
opti.set_initial(U, zeros(1,N_phase));

% restrict velocity


% for solution with guard

%load('relaxed_solution_controlled_4.mat')
% load('relaxed_solution_2.mat')
% load('strict_solution_5.mat')
% load('kinda_strict_solution_controlled_1.mat')
% opti.set_initial(X, x_sol);
% opti.set_initial(T_fl, time_fl);
% opti.set_initial(T_st, time_st);
% opti.set_initial(U, zeros(1,N_phase));

% show progress of optimization 
opti.callback(@(i) plot(opti.debug.value(X(1,1:N+1)), opti.debug.value(X(2,1:N+1)),'b') )
grid on;

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
time_st_L1 = sol.value(T_st_L1)
time_st_L2 = sol.value(T_st_L2)
x_sol = sol.value(X)
u_sol = sol.value(U)
phi0_sol = sol.value(phi0)
theta0_sol = sol.value(theta0)
alpha0_sol = sol.value(alpha0)
omega_sol = sol.value(omega)


%% Plot 

x_hip_sol = sol.value(X(1,:));
y_hip_sol = sol.value(X(2,:));
phi_sol = sol.value(X(3,:));
phi1_sol = sol.value(X(4,:));
theta_sol = sol.value(X(5,:));
theta1_sol = sol.value(X(6,:));
angle1_sol = phi_sol - phi1_sol;%%%%%%%%%%
angle2_sol = theta_sol - theta1_sol;%%%%%%%%%%
phi_td = sol.value(phi_sol(N_phase+1)); 
theta_td = sol.value(theta_sol(N_phase+1)); 
angle1_td = sol.value(phi - phi1);%%%%%%%%%%%
angle2_td = sol.value(theta - theta1);%%%%%%%%%%%
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
