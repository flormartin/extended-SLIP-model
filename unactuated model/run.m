clc; clear all; %close all;

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
c1 = 20 * 1000;


c1 = 20*1000;
alpha0 = (90 - 62) * pi/ 180.;
phi0 = -35*pi/180;
omega = 50*pi/180; %rad/s

x0 = [0.; .95; -35*pi/180; 2; 2; 0];
t_apex = x0(5)/g;

opts_flight = odeset('Events', @guard_flight);
opts_stance = odeset('Events', @guard_stance);

step_width = 0.001;
tend = 2.;
tspan = 0:step_width:tend;

x_full = [];
t_full = [];

% one flight phase and one stance phase
[tout, xout, te, xe, ie] = ode45(@(t,x)flight(t,x,alpha0,omega), tspan, x0, opts_flight);

t_full = [t_full; tout + 0];
x_full = [x_full; xout, zeros(size(xout,1),1)];
y_full = [t_full; tout];

offset = xe(1) + L0*sin(xe(3));
figure;plot(offset,0,'x','Color','r','Markersize',15,'linewidth',2); hold on
x0 = xe';
tspan = te:step_width:tend;

[tout, xout, te, xe, ie] = ode45(@(t,x)stance(x,phi0,offset), tspan, x0, opts_stance);

t_full = [t_full; tout + t_full(end)];
theta = xout(:,3) - atan((offset-(xout(:,1)+sin(xout(:,3))*l1))./(xout(:,2)-cos(xout(:,3))*l1));
x_full = [x_full; xout,theta];
y_full = [t_full; tout];

% x0 = xe';
% tspan = te:step_width:tend;
% t_apex = te + xe(5)/g;
% 
% [tout, xout, te, xe, ie] = ode45(@(t,x)flight(t,x,alpha0,omega), tspan, x0, opts_flight);
% 
% t_full = [t_full; tout + 0];
% x_full = [x_full; xout, zeros(size(xout,1),1)];
% y_full = [t_full; tout];
% 
% offset = xe(1) + L0*sin(xe(3));
% plot(offset,0,'x','Color','r','Markersize',15,'linewidth',2); hold on
% x0 = xe';
% tspan = te:step_width:tend;
% 
% [tout, xout, te, xe, ie] = ode45(@(t,x)stance(x,phi0,offset), tspan, x0, opts_stance);
% 
% t_full = [t_full; tout + t_full(end)];
% theta = xout(:,3) - atan((offset-(xout(:,1)+sin(xout(:,3))*l1))./(xout(:,2)-cos(xout(:,3))*l1));
% x_full = [x_full; xout,theta];
% y_full = [t_full; tout];
% 
% x0 = xe';
% tspan = te:step_width:tend;
% t_apex = te + xe(5)/g;

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
figure;plot(x_full(:,7)*180/pi)
