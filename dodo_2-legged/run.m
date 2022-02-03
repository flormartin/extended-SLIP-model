clc; clear all; close all;

global L0 l0 l1 M m_B m_L alpha0 g c_hip dhip c_knee phi1_0 theta1_0 dknee c1 offset t_apex omega phi0 theta0;
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
phi0 = 0;
theta0 = 0;
phi1_0 = 0;
theta1_0 = 0;

c1 = 3000;
alpha0 = (90 - 62) * pi/ 180.;
t_apex = 0;
omega = 14*pi/180; %rad/s
x0 = [0.; 0.4; alpha0; 0; -alpha0; 0; 2; -0.2; 0; 0; 0; 0];


opts_stance_L1 = odeset('Events', @guard_stance_leg2);
opts_stance_L2 = odeset('Events', @guard_stance_leg1);

step_width = 0.001;
tend = 10.;
tspan = 0:step_width:tend;

x_full = [];
t_full = [];

%% one stance foot 1 and one stance foot 2

%Stance-phase foot 1, phi-angle
offset = x0(1) + L0*sin(x0(3));
[tout, xout, te, xe, ie] = ode45(@(t,x)mode1_2L_run(t,x), tspan, x0, opts_stance_L1);
t_full = [t_full; tout + 0];
x_full = [x_full; xout];
y_full = [t_full; tout];

plot(offset,0,'x','Color','r','Markersize',15,'linewidth',2); hold on

%Stance-phase foot 2, theta-angle
x0 = xe';
tspan = te:step_width:tend;
t_apex = (te - t_apex) / 2;
offset = xe(1) + L0*sin(xe(5));

[tout, xout, te, xe, ie] = ode45(@(t,x)mode2_2L_run(t,x), tspan, x0, opts_stance_L2);
t_full = [t_full; tout + t_full(end)];
x_full = [x_full; xout];
y_full = [t_full; tout];

x0 = xe';
tspan = te:step_width:tend;

%% plots

for i = 1:length(t_full)
    if mod(i,20) == 0
        plot(x_full(i,1), x_full(i,2),'o','Color','k')%hip
        plot(x_full(i,1)+l1*sin(x_full(i,3)), x_full(i,2)-l1*cos(x_full(i,3)),'o','Color','k')%knee of phi-foot
        plot(x_full(i,1)+l1*sin(x_full(i,5)), x_full(i,2)-l1*cos(x_full(i,5)),'o','Color','blue')%knee of theta-foot
        plot(x_full(i,1)+l1*sin(x_full(i,3))+l0*sin(x_full(i,3)-x_full(i,4)), x_full(i,2)-l1*cos(x_full(i,3))-l0*cos(x_full(i,3)-x_full(i,4)),'o','Color','k')%foot of phi
        plot(x_full(i,1)+l1*sin(x_full(i,5))+l0*sin(x_full(i,5)-x_full(i,6)), x_full(i,2)-l1*cos(x_full(i,5))-l0*cos(x_full(i,5)-x_full(i,6)),'o','Color','blue')%foot of theta
        plot([x_full(i,1),x_full(i,1)+l1*sin(x_full(i,3))],[x_full(i,2),x_full(i,2)-l1*cos(x_full(i,3))],'Color','k')
        plot([x_full(i,1),x_full(i,1)+l1*sin(x_full(i,5))],[x_full(i,2),x_full(i,2)-l1*cos(x_full(i,5))],'Color','blue')
        plot([x_full(i,1)+l1*sin(x_full(i,3)),x_full(i,1)+l1*sin(x_full(i,3))+l0*sin(x_full(i,3)-x_full(i,4))],[x_full(i,2)-l1*cos(x_full(i,3)),x_full(i,2)-l1*cos(x_full(i,3))-l0*cos(x_full(i,3)-x_full(i,4))],'Color','k')
        plot([x_full(i,1)+l1*sin(x_full(i,5)),x_full(i,1)+l1*sin(x_full(i,5))+l0*sin(x_full(i,5)-x_full(i,6))],[x_full(i,2)-l1*cos(x_full(i,5)),x_full(i,2)-l1*cos(x_full(i,5))-l0*cos(x_full(i,5)-x_full(i,6))],'Color','blue')
    end
end
figure;plot(x_full(:,1),x_full(:,2))
title('hip-position')
