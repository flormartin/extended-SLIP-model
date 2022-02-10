clc; clear all; close all;

%parameters
global L0 l0 l1 M m_B m_L g c_phi d1 c_theta d2 c1 offset t_apex;
L0=1; M=80; g=9.81;
l0 = .6*L0;                     %lower leg length
l1 = .4*L0;                     %upper leg length
m_L = .32*M;                    %leg mass
m_B = M-m_L;                  %rest of body mass
c_phi = 750;                    %Nm/rad
d1 = 2*sqrt(c_phi*m_L);
c_theta = 1000;
d2 = 2*sqrt(c_phi*m_L);

load('almoststrict.mat')

c1 = 20*1000;
alpha0 = alpha0_sol;
phi0 = phi0_sol;
omega = omega_sol; %rad/s
x0 = x_sol(:,1);
% x0(4)=x0(4)+.9;

% guard functions
opts_flight = odeset('Events', @guard_flight);
opts_stance = odeset('Events', @guard_stance);


t_apex = x0(5)/g;
step_width = 0.001;
tend = 2.;
tspan = 0:step_width:tend;

x_full = [];
t_full = [];
steps = 0;

%flight phase
[tout, xout, te, xe, ie] = ode45(@(t,x)mode2(t,x,alpha0,omega), tspan, x0, opts_flight);

t_full = [t_full; tout + 0];
x_full = [x_full; xout, zeros(size(xout,1),1)];
y_full = [t_full; tout];

offset = xe(1) + L0*sin(xe(3));
figure;plot(offset,0,'x','Color','r','Markersize',15,'linewidth',2); hold on
x0 = xe';
tspan = tout(end):step_width:tend;

while(isempty(xe)==false)
%stance phase
[tout, xout, te, xe, ie] = ode45(@(t,x)mode1(x,phi0,offset), tspan, x0, opts_stance);

t_full = [t_full; tout + t_full(end)];
theta = xout(:,3) - atan((offset-(xout(:,1)+sin(xout(:,3))*l1))./(xout(:,2)-cos(xout(:,3))*l1));
x_full = [x_full; xout,theta];
y_full = [t_full; tout];

if isempty(xe)
    break
end
x0 = xout(end,:)';
tspan = tout(end):step_width:tend;
t_apex = x0(5)/g;

%flight phase
[tout, xout, te, xe, ie] = ode45(@(t,x)mode2(t,x,alpha0,omega), tspan, x0, opts_flight);

t_full = [t_full; tout + 0];
x_full = [x_full; xout, zeros(size(xout,1),1)];
y_full = [t_full; tout];

if isempty(xe)
    break
end
offset = xout(end,1) + L0*sin(xout(end,3));
plot(offset,0,'x','Color','r','Markersize',15,'linewidth',2); hold on
x0 = xout(end,:)';
tspan = tout(end):step_width:tend;
steps = steps + 1;
end


%% plots
for i = 1:length(t_full)
    if mod(i,5) == 0
        plot(x_full(i,1), x_full(i,2),'o','Color','k')%hip
        plot(x_full(i,1)+l1*sin(x_full(i,3)), x_full(i,2)-l1*cos(x_full(i,3)),'o','Color','k')%knee
        plot(x_full(i,1)+l1*sin(x_full(i,3))+l0*sin(x_full(i,3)-x_full(i,7)), x_full(i,2)-l1*cos(x_full(i,3))-l0*cos(x_full(i,3)-x_full(i,7)),'o','Color','k')%foot
        plot([x_full(i,1),x_full(i,1)+l1*sin(x_full(i,3))],[x_full(i,2),x_full(i,2)-l1*cos(x_full(i,3))],'Color','k')
        plot([x_full(i,1)+l1*sin(x_full(i,3)),x_full(i,1)+l1*sin(x_full(i,3))+l0*sin(x_full(i,3)-x_full(i,7))],[x_full(i,2)-l1*cos(x_full(i,3)),x_full(i,2)-l1*cos(x_full(i,3))-l0*cos(x_full(i,3)-x_full(i,7))],'Color','k')
    end
end

figure;plot(x_full(:,2),x_full(:,5))

% figure;plot(x_full(:,7)*180/pi)
