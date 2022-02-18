clc; clear all; close all;

%load solution
load('strict3.mat')

[ap,i] = min(abs(x_sol(5,41:end)));
i = i + 40;
x0 = [0;x_sol(2:end,i)]';

%limit cycle
figure; hold on
plot(x_sol(2,:),x_sol(5,:),'color','#0072BD','linewidth',2)
plot(x_sol(2,i),x_sol(5,i),'.','color','k','markersize',20)
figure; hold on
plot(x_sol(3,:),x_sol(6,:),'color','#0072BD','linewidth',2)
plot(x_sol(3,i),x_sol(6,i),'.','color','k','markersize',20)

%calculate dP and eigenvalues

% poincare(x0,phi0_sol,alpha0_sol, omega_sol)
% x0

epsilon = 1e-5;
dP = [];

for i = 2:6
    r = zeros(1,6);
    r(i) = epsilon;
    df = (poincare(x0+r,phi0_sol,alpha0_sol,omega_sol) - poincare(x0-r,phi0_sol,alpha0_sol,omega_sol))...
        /(2.*epsilon);
    dP(i-1,:) = df(2:6);
end

dP = dP;
eigenvalues = eig(dP)


function [P, x_full, te] = poincare(x0, phi0, alpha0, omega)

% Parameters
global L0 l0 l1 M m_B m_L g c_phi d1 c_theta d2 c1 t_apex offset;
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

opts_flight = odeset('Events', @guard_flight);
opts_stance = odeset('Events', @guard_stance);
opts_apex = odeset('Events', @guard_apex);
step_width=.001;
x_full = [];
t_full = [];

t_apex = x0(5)/g;
tspan = 0:step_width:4.0;
[tout, xout, te, xe, ie] = ode45(@(t,x)mode2(t,x,alpha0,omega), tspan, x0, opts_flight);
x_full = [x_full; xout];
t_full = [t_full; tout];

tspan = te:step_width:4.0;
offset = xe(1) + L0 * sin(xe(3));
x0 = xe';
[tout, xout, te, xe, ie] = ode45(@(t,x)mode1st(t,x,phi0,offset), tspan, x0, opts_stance);
x_full = [x_full; xout];
t_full = [t_full; tout];

t_apex = te+xe(5)/g;
tspan = te:step_width:4.0;
x0 = xe';
[tout, xout, te, xe, ie] = ode45(@(t,x)mode2(t,x,alpha0,omega), tspan, x0, opts_apex);
x_full = [x_full; xout];
t_full = [t_full; tout];

% figure;plot(x_full(:,1),x_full(:,2))
% figure;plot(x_full(:,2),x_full(:,5))
% figure;plot(x_full(:,3),x_full(:,6))

P = xe;
end