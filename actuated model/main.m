clc; clear all; close all;
addpath('C:\Program Files\MATLAB\R2021a\toolbox\casadi-windows-matlabR2016a-v355')
import casadi.*

%% Initialization

opti = casadi.Opti();

% decision_variables = [X, T_flight, T_stance]

% Number of nodes
N = 80;
N_phase = N/2;

% x = [x y phi dx dy dphi]
X = opti.variable(6, N+1);
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
m_L = .32*M;                    %leg mass
m_B = M-m_L;                    %rest of body mass
c_phi = 750;                    %Nm/rad
d1 = 2*sqrt(c_phi*m_L);
c_theta = 1000;
d2 = 2*sqrt(c_theta*m_L);
c1 = 20 * 1000;

c1 = 20*1000;
%alpha0 = (90 - 62) * pi/ 180.;
t_apex = X(5,1)/g;              %vertical initial speed divided by g
%omega = 50*pi/180; %rad/s

%% Flight phase - Mode 2

h_fl = T_fl/N_phase;    % step width

% Runge-Kutta
for k = 1:N_phase
    
    t = k * h_fl;
    k1 = mode2(t, X(:,k), alpha0, omega);
    k2 = mode2(t, X(:,k) + 0.5 * h_fl * k1, alpha0, omega);
    k3 = mode2(t, X(:,k) + 0.5 * h_fl * k2, alpha0, omega);
    k4 = mode2(t, X(:,k) + h_fl * k3, alpha0, omega);
    
    x_next = X(:,k) + 1./6. * h_fl * (k1 + 2 * k2 + 2 * k3 + k4);
    
    opti.subject_to(X(:,k+1) == x_next);
end

phi = X(3,N_phase+1);
offset = X(1,N_phase+1) + L0 * sin(phi);

%% Stance phase - Mode 1

h_st = T_st/N_phase;    % step width

% Runge-Kutta
for k = N_phase+1:N

    k1 = mode1(X(:,k), phi0, offset);
    k2 = mode1(X(:,k) + 0.5 * h_st * k1, phi0, offset);
    k3 = mode1(X(:,k) + 0.5 * h_st * k2, phi0, offset);
    k4 = mode1(X(:,k) + h_st * k3, phi0, offset);
    
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

%% Guard conditions

% liftoff
x_lo_knee = X(1,N+1) - offset + l1 * sin(X(3,N+1)); 
y_lo_knee = X(2,N+1) - l1 * cos(X(3,N+1)) ;
% constraint: length of lower leg equal to 0.6
opti.subject_to(sqrt(x_lo_knee^2 + y_lo_knee^2) == l0);  

% touchdown: height of hip = cos(phi)*L0 
opti.subject_to(X(2,N_phase+1) == L0 * cos(X(3,N_phase+1)));

%% Path contraints
% x = [x y phi theta dx dy dphi dtheta]
opti.subject_to(T_fl >= 0.2);  
opti.subject_to(T_st >= 0.2);

opti.subject_to(X(1,1) == 0); % starting point - horizontal coordinate x
opti.subject_to(X(2,:) >= .4); % height of hip
opti.subject_to(X(3,:) >= -pi/2);   % phi 
opti.subject_to(X(3,:) <= pi/2);
opti.subject_to(X(3,1) <= 0);
% opti.subject_to(X(3,N+1) == atan((offset-(X(1,N+1)+sin(X(3,N+1))*l1))/(X(2,N+1)-cos(X(3,N+1))*l1))); %leads to very small step solution
opti.subject_to(X(4,:) >= 0);   % dx - hip keeps moving forward not backward
opti.subject_to(X(5,1) > 0);    % hip moving upwards at liftoff

% new constraints including the actuation parameters as optimization
% variables
opti.subject_to(phi0 <= 0);
opti.subject_to(phi0 >= -pi/2);
opti.subject_to(alpha0 <= pi/2);
opti.subject_to(alpha0 >= -pi/2);
opti.subject_to(omega >= 0);

%% Initial guess

% for relaxed solution - without guard 
% x0 = [0.; 1.05; -35 * pi/ 180.; 2.5; 2; 0];
% % x0 = ones(6,1)*.5;
% opti.set_initial(X, repmat(x0,1,N+1));
% opti.set_initial(phi0, -35*pi/180);
% opti.set_initial(alpha0, -35*pi/180);
% opti.set_initial(omega, 50*pi/180);
% opti.set_initial(T_fl, 0.2);
% opti.set_initial(T_st, 0.2);

% for solution with guard

load('relaxed.mat')
opti.set_initial(X, x_sol);
opti.set_initial(T_fl, T_fl_sol);
opti.set_initial(T_st, T_st_sol);
opti.set_initial(phi0, phi0_sol);
opti.set_initial(alpha0, alpha0_sol);
opti.set_initial(omega, omega_sol);

% show progress of optimization 
opti.callback(@(i) plot(opti.debug.value(X(1,:)), opti.debug.value(X(2,:)),'b') )

options = struct;                                                     % alternative 2
options.print_time = false;
options.ipopt.max_iter = 800;
opti.solver('ipopt',options);   % interior point method
sol = opti.solve();

%% Plot 

x_hip_sol = sol.value(X(1,:));
y_hip_sol = sol.value(X(2,:));
phi_sol = sol.value(X(3,:));
offset_sol = x_hip_sol(N_phase+1) + L0 * sin(phi_sol(N_phase+1));



% Knee
x_knee = x_hip_sol + l1 * sin(phi_sol);
y_knee = y_hip_sol - l1 * cos(phi_sol);

% Foot
l_s_sol(1:N_phase) = l0;
l_s_sol(N_phase+1:N+1) = sqrt((offset_sol-(x_hip_sol(N_phase+1:N+1)+l1*sin(phi_sol(N_phase+1:N+1)))).^2+(y_hip_sol(N_phase+1:N+1)-l1*cos(phi_sol(N_phase+1:N+1))).^2);

angle_sol = phi_sol(1:N_phase);
angle_sol(N_phase+1:N+1)=atan((offset_sol-(x_knee(N_phase+1:N+1)))./y_knee(N_phase+1:N+1));
theta_sol = zeros(1,N+1);
theta_sol(N_phase+1:N+1)=phi_sol(N_phase+1:N+1)-angle_sol(N_phase+1:N+1);

x_foot = x_knee + l_s_sol .* sin(angle_sol);
y_foot = y_knee - l_s_sol .* cos(angle_sol);

%% 
figure; hold on;
for i = 1:N+1
    if mod(i,1)==0
        plot(x_hip_sol(i), y_hip_sol(i),'o','Color','k');   % position of hip
        plot(x_knee(i), y_knee(i),'o','Color','b'); % position of knee
        plot(x_foot(i), y_foot(i),'o','Color','k');  % position of foot
        plot([x_hip_sol(i),x_knee(i)],[y_hip_sol(i),y_knee(i)],'Color','k');    % upper leg
        plot([x_knee(i),x_foot(i)],[y_knee(i),y_foot(i)],'Color','k');  % lower leg
    end
end
plot(x_foot, y_foot,'-','Color','r');  %  foot
plot(x_knee, y_knee,'-','Color','r');  %  knee
plot(x_hip_sol, y_hip_sol,'-','Color','r');  %  hip

%% Save solution

x_sol = sol.value(X);
T_fl_sol = sol.value(T_fl);
T_st_sol = sol.value(T_st);
phi0_sol = sol.value(phi0);
alpha0_sol = sol.value(alpha0);
omega_sol = sol.value(omega);
% 
% % x_sol = opti.debug.value(X);
% % T_fl_sol = opti.debug.value(T_fl);
% % T_st_sol = opti.debug.value(T_st);
% % phi0_sol = opti.debug.value(phi0);
% % alpha0_sol = opti.debug.value(alpha0);
% % omega_sol = opti.debug.value(omega);