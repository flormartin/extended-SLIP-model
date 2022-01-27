
%% Initialization
% 
% % alpha0_ = sol.value(alpha0);
% global offset l0 l1;
% L0=1; M=80; g=9.81;
% l0 = .6*L0;                     %lower leg length
% l1 = .4*L0;                     %upper leg length
load('Nice_Solution_01.mat');
% Number of nodes
N = 80;
N_phase = N/2;
alpha0 = alpha0_sol;
omega = omega_sol;
phi0 = phi0_sol;


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
%alpha0 = (90 - 62) * pi/ 180.;
t_apex = x_sol(6,1);
%omega = 50*pi/180; %rad/s

X = x_sol;


%% Flight phase - Mode 2

h_fl_ = time_fl/N_phase;    % step width

% Runge-Kutta
for k = 1:N_phase
    
    t = k * h_fl_;
    k1 = mode2(t, X(:,k), alpha0, omega);
    k2 = mode2(t, X(:,k) + 0.5 * h_fl_ * k1, alpha0, omega);
    k3 = mode2(t, X(:,k) + 0.5 * h_fl_ * k2, alpha0, omega);
    k4 = mode2(t, X(:,k) + h_fl_ * k3, alpha0, omega);
    
    X(:,k+1) = X(:,k) + 1./6. * h_fl_ * (k1 + 2 * k2 + 2 * k3 + k4);
    
end

phi = X(3,N_phase+1);
theta = X(4,N_phase+1);
offset = X(1,N_phase+1) + l1 * sin(phi) + l0 * sin(phi - theta);

%% Stance phase - Mode 1

h_st_ = time_st/N_phase;    % step width


% Runge-Kutta
for k = N_phase+1:N
    j = k-N_phase;

    k1 = mode1(X(:,k), phi0, offset, u_sol(j));
    k2 = mode1(X(:,k) + 0.5 * h_st_ * k1, phi0, offset, u_sol(j));
    k3 = mode1(X(:,k) + 0.5 * h_st_ * k2, phi0, offset, u_sol(j));
    k4 = mode1(X(:,k) + h_st_ * k3, phi0, offset, u_sol(j));
    
    X(:,k+1) = X(:,k) + 1./6. * h_st_ * (k1 + 2 * k2 + 2 * k3 + k4);
    
end

%% next round
X2(:,1) = X(:,N+1)

%% Flight phase - Mode 2

h_fl_ = time_fl/N_phase;    % step width

% Runge-Kutta
for k = 1:N_phase
    
    t = k * h_fl_;
    k1 = mode2(t, X2(:,k), alpha0, omega);
    k2 = mode2(t, X2(:,k) + 0.5 * h_fl_ * k1, alpha0, omega);
    k3 = mode2(t, X2(:,k) + 0.5 * h_fl_ * k2, alpha0, omega);
    k4 = mode2(t, X2(:,k) + h_fl_ * k3, alpha0, omega);
    
    X2(:,k+1) = X2(:,k) + 1./6. * h_fl_ * (k1 + 2 * k2 + 2 * k3 + k4);
    
end

phi = X2(3,N_phase+1);
theta = X2(4,N_phase+1);
offset = X2(1,N_phase+1) + l1 * sin(phi) + l0 * sin(phi - theta);

%% Stance phase - Mode 1

h_st_ = time_st/N_phase;    % step width


% Runge-Kutta
for k = N_phase+1:N
    j = k-N_phase;

    k1 = mode1(X2(:,k), phi0, offset, u_sol(j));
    k2 = mode1(X2(:,k) + 0.5 * h_st_ * k1, phi0, offset, u_sol(j));
    k3 = mode1(X2(:,k) + 0.5 * h_st_ * k2, phi0, offset, u_sol(j));
    k4 = mode1(X2(:,k) + h_st_ * k3, phi0, offset, u_sol(j));
    
    X2(:,k+1) = X2(:,k) + 1./6. * h_st_ * (k1 + 2 * k2 + 2 * k3 + k4);
    
end

%% next round
X3(:,1) = X2(:,N+1)

%% Flight phase - Mode 2

h_fl_ = time_fl/N_phase;    % step width

% Runge-Kutta
for k = 1:N_phase
    
    t = k * h_fl_;
    k1 = mode2(t, X3(:,k), alpha0, omega);
    k2 = mode2(t, X3(:,k) + 0.5 * h_fl_ * k1, alpha0, omega);
    k3 = mode2(t, X3(:,k) + 0.5 * h_fl_ * k2, alpha0, omega);
    k4 = mode2(t, X3(:,k) + h_fl_ * k3, alpha0, omega);
    
    X3(:,k+1) = X3(:,k) + 1./6. * h_fl_ * (k1 + 2 * k2 + 2 * k3 + k4);
    
end

phi = X3(3,N_phase+1);
theta = X3(4,N_phase+1);
offset = X3(1,N_phase+1) + l1 * sin(phi) + l0 * sin(phi - theta);

%% Stance phase - Mode 1

h_st_ = time_st/N_phase;    % step width


% Runge-Kutta
for k = N_phase+1:N
    j = k-N_phase;

    k1 = mode1(X3(:,k), phi0, offset, u_sol(j));
    k2 = mode1(X3(:,k) + 0.5 * h_st_ * k1, phi0, offset, u_sol(j));
    k3 = mode1(X3(:,k) + 0.5 * h_st_ * k2, phi0, offset, u_sol(j));
    k4 = mode1(X3(:,k) + h_st_ * k3, phi0, offset, u_sol(j));
    
    X3(:,k+1) = X3(:,k) + 1./6. * h_st_ * (k1 + 2 * k2 + 2 * k3 + k4);
    
end

%% Plot

figure;
plot(X(1,:),X(2,:))
hold on;
plot(X2(1,:),X2(2,:))
plot(X3(1,:),X3(2,:))
X_ = [X(:,1:end-1), X2(:,1:end-1), X3(:,1:end)];

%% Plot 

x_hip_sol = X_(1,:);
y_hip_sol = X_(2,:);
phi_sol = X_(3,:);
theta_sol = X_(4,:);
angle_sol = phi_sol - theta_sol;%%%%%%%%%%
phi_td = phi_sol(N_phase+1); 
theta_td = theta_sol(N_phase+1); 
angle_td = phi - theta;%%%%%%%%%%%
offset_sol = x_hip_sol(N_phase+1) + l_u * sin(phi_td) + l_l * sin(angle_td);

% Knee
x_knee = zeros(1,3*N+1);
y_knee = zeros(1,3*N+1);
x_knee = x_hip_sol + l_u * sin(phi_sol);
y_knee = y_hip_sol - l_u * cos(phi_sol);

% Foot 
x_foot = zeros(1,3*N+1);
y_foot = zeros(1,3*N+1);
l_s_sol = zeros(1,3*N+1);
l_s_sol(1:N_phase) = 0.6;
l_s_sol(N+1:N+N_phase) = 0.6;
l_s_sol(N+N+1:N+N+N_phase) = 0.6;
offset_hip = X(1,end);
for i = N_phase+1:N+1
    l_s_sol(i) = sqrt((offset_sol-(x_hip_sol(i)+l1*sin(phi_sol(i))))^2+(y_hip_sol(i)-l1*cos(phi_sol(i)))^2);
end
for i = N+N_phase+1:N+N+1
    l_s_sol(i) = sqrt((offset_hip + offset_sol-(x_hip_sol(i)+l1*sin(phi_sol(i))))^2+(y_hip_sol(i)-l1*cos(phi_sol(i)))^2);
end
for i = N+N+N_phase+1:N+N+N+1
    l_s_sol(i) = sqrt((2 * offset_hip + offset_sol-(x_hip_sol(i)+l1*sin(phi_sol(i))))^2+(y_hip_sol(i)-l1*cos(phi_sol(i)))^2);
end

x_foot = x_knee + l_s_sol .* sin(angle_sol);
y_foot = y_knee - l_s_sol .* cos(angle_sol);

 
figure; hold on;
plot(x_foot, y_foot,'-','Color','r');  %  foot
plot(x_knee, y_knee,'-','Color','r');  %  knee
plot(x_hip_sol, y_hip_sol,'-','Color','r');  %  hip
for i = 1:3*N+1
    plot(x_hip_sol(i), y_hip_sol(i),'o','Color','k');   % position of hip
    plot(x_knee(i), y_knee(i),'o','Color','b'); % position of knee
    plot(x_foot(i), y_foot(i),'o','Color','k');  % position of foot
    plot([x_hip_sol(i),x_knee(i)],[y_hip_sol(i),y_knee(i)],'Color','k');    % upper leg
    plot([x_knee(i),x_foot(i)],[y_knee(i),y_foot(i)],'Color','k');  % lower leg
end

