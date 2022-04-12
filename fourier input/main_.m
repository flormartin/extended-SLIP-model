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
n_coeff = 1;
a = opti.variable(n_coeff+1);
b = opti.variable(n_coeff);

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

c1 = 20*1000;
%alpha0 = (90 - 62) * pi/ 180.;

%omega = 50*pi/180; %rad/s

%% Stance phase - Mode 1

offset = X(1,1) + L0 * sin(X(3,1));

h_st = T_st/N_phase;    % step width

% Runge-Kutta
for k = 1:N_phase
    
    t = k*h_st;
    i = (1:n_coeff)';
    phi0 = a(1) + sum(a(i+1).*cos(i*t*2*pi/(T_st+T_fl))+b(i).*sin(i*t*2*pi/(T_st+T_fl)));
    
    k1 = mode1(X(:,k), phi0, offset);
    k2 = mode1(X(:,k) + 0.5 * h_st * k1, phi0, offset);
    k3 = mode1(X(:,k) + 0.5 * h_st * k2, phi0, offset);
    k4 = mode1(X(:,k) + h_st * k3, phi0, offset);
    
    x_next = X(:,k) + 1./6. * h_st * (k1 + 2 * k2 + 2 * k3 + k4);
    
    opti.subject_to(X(:,k+1) == x_next);
end

%% Flight phase - Mode 2

t_apex = T_st + X(5,1)/g;              %vertical initial speed divided by g

h_fl = T_fl/N_phase;    % step width

% Runge-Kutta
for k = N_phase+1:N
    
    t = T_st + (k-N_phase) * h_fl;
    i = (1:n_coeff)';
    phi0 = a(1) + sum(a(i+1).*cos(i*t*2*pi/(T_st+T_fl))+b(i).*sin(i*t*2*pi/(T_st+T_fl)));
    
    k1 = mode2(X(:,k), phi0);
    k2 = mode2(X(:,k) + 0.5 * h_fl * k1, phi0);
    k3 = mode2(X(:,k) + 0.5 * h_fl * k2, phi0);
    k4 = mode2(X(:,k) + h_fl * k3, phi0);
    
    x_next = X(:,k) + 1./6. * h_fl * (k1 + 2 * k2 + 2 * k3 + k4);
    
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
x_lo_knee = X(1,N_phase+1) - offset + l1 * sin(X(3,N_phase+1)); 
y_lo_knee = X(2,N_phase+1) - l1 * cos(X(3,N_phase+1)) ;
% constraint: length of lower leg equal to 0.6
opti.subject_to(x_lo_knee^2 + y_lo_knee^2 == l0^2);  
% 
% touchdown: height of hip = cos(phi)*L0 
opti.subject_to(X(2,1) == L0 * cos(X(3,1)));

%% Path contraints
% x = [x y phi theta dx dy dphi dtheta]
opti.subject_to(T_fl >= 0.15);  
opti.subject_to(T_st >= 0.15);

opti.subject_to(X(1,1) == 0);           % starting point - horizontal coordinate x
opti.subject_to(X(2,1) >= 0.9); 
opti.subject_to(0.1 < X(1,N+1) < 5);    % length of step
opti.subject_to(0.4 < X(2,:) < 2);      % height of hip
opti.subject_to(-pi/2 < X(3,:) < pi/2); % phi 
opti.subject_to(X(4,:) >= 0);           % dx - hip keeps moving forward not backward
% opti.subject_to(X(4,1) == 5);           % initial speed
opti.subject_to(3.5 < X(4,1) < 6);      % initial speed range
% opti.subject_to(mean(X(4,:)) == 3.6149);
% opti.subject_to(abs(X(6,:)) <= 100); %limit rotation speed around hip

%theta constraints
% opti.subject_to(X(3,1)-atan((offset-(X(1,1)+sin(X(3,1))*l1))/(X(2,1)-cos(X(3,1))*l1))==0);
opti.subject_to(X(3,:)-atan((offset-(X(1,:)+sin(X(3,:))*l1))./(X(2,:)-cos(X(3,:))*l1))>=0);
% opti.subject_to(X(3,N_phase+1)-atan((offset-(X(1,N_phase+1)+sin(X(3,N_phase+1))*l1))/(X(2,N_phase+1)-cos(X(3,N_phase+1))*l1))==0);

%% Initial guess

% % for relaxed solution - without guard 
% x0 = [0.; .8; 68 * pi/ 180.; 5; -1; 0];
% opti.set_initial(X, repmat(x0,1,N+1));
% opti.set_initial(a,ones(1,n_coeff+1))
% opti.set_initial(b,ones(1,n_coeff))
% opti.set_initial(T_fl, 0.2);
% opti.set_initial(T_st, 0.2);

% for solution with guard
load('relaxed2_.mat')
opti.set_initial(X, x_sol);
opti.set_initial(T_fl, T_fl_sol);
opti.set_initial(T_st, T_st_sol);
opti.set_initial(a, a_sol);
opti.set_initial(b, b_sol);

% show progress of optimization 
opti.callback(@(i) plot(opti.debug.value(X(1,:)), opti.debug.value(X(2,:)),'b') )

options = struct;
options.print_time = false;
options.ipopt.max_iter = 4000;
opti.solver('ipopt',options);   % interior point method
sol = opti.solve();

%% Plot Solution

x_hip_sol = sol.value(X(1,:));
y_hip_sol = sol.value(X(2,:));
phi_sol = sol.value(X(3,:));
offset_sol = x_hip_sol(1) + L0 * sin(phi_sol(1));

% Knee
x_knee = x_hip_sol + l1 * sin(phi_sol);
y_knee = y_hip_sol - l1 * cos(phi_sol);

% Foot
l_s_sol(1:N_phase) = sqrt((offset_sol-(x_hip_sol(1:N_phase)+l1*sin(phi_sol(1:N_phase)))).^2+(y_hip_sol(1:N_phase)-l1*cos(phi_sol(1:N_phase))).^2);
l_s_sol(N_phase+1:N+1) = l0;

% % no retraction of lower leg
% angle_sol(1:N_phase)=atan((offset_sol-(x_knee(1:N_phase)))./y_knee(1:N_phase));
% angle_sol(N_phase+1:N+1) = phi_sol(N_phase+1:N+1);

% linear retraction of lower leg
angle_sol(1:N_phase)=atan((offset_sol-(x_knee(1:N_phase)))./y_knee(1:N_phase));
angle_sol(N_phase+1:N+1) = linspace(angle_sol(N_phase),phi_sol(N+1),N_phase+1);

theta_sol = zeros(1,N+1);
theta_sol(1:N_phase)=phi_sol(1:N_phase)-angle_sol(1:N_phase);

x_foot = x_knee + l_s_sol .* sin(angle_sol);
y_foot = y_knee - l_s_sol .* cos(angle_sol);

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

k = 1:n_coeff;
T = sol.value(T_fl+T_st);
t_step(1:N_phase) = (sol.value(h_st):sol.value(h_st):sol.value(T_st));
t_step(N_phase+1:N+1) = (sol.value(T_st):sol.value(h_fl):T);
fourier_phi0(1:N+1) = sol.value(a(1))*ones(size(t_step)) + sum(sol.value(a(k+1))*ones(size(t_step)).*cos(k'*t_step*2*pi/T)+sol.value(b(k))*ones(size(t_step)).*sin(k'*t_step*2*pi/T),1);
fourier_phi0(N+2) = fourier_phi0(1);
t_step(N+2) = t_step(end)+t_step(1);

figure;
plot(t_step,fourier_phi0)

%% Save solution

x_sol = sol.value(X);
T_fl_sol = sol.value(T_fl);
T_st_sol = sol.value(T_st);
a_sol = sol.value(a);
b_sol = sol.value(b);

% % x_sol = opti.debug.value(X);
% % T_fl_sol = opti.debug.value(T_fl);
% % T_st_sol = opti.debug.value(T_st);
% % phi0_sol = opti.debug.value(phi0);
% % alpha0_sol = opti.debug.value(alpha0);
% % omega_sol = opti.debug.value(omega);
