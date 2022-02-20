clear all; close all;

%load solution
load('strict1902_1200_40k_7500.mat')
a = a_sol;
b = b_sol;
T = T_st_sol+T_fl_sol;

[ap,i] = min(abs(x_sol(5,41:end)));
x0 = [0;x_sol(2:end,i + 40)]';
toffset = T_st_sol + i * T_fl_sol/41;
steps = 0;

% Parameters
global L0 l0 l1 M m_B m_L g c_phi d1 c_theta d2 c1 t_apex offset;
L0=1; M=80; g=9.81;
l0 = .6*L0;                     %lower leg length
l1 = .4*L0;                     %upper leg length
m_L = .32*M;                    %leg mass
m_B = M-m_L;                    %rest of body mass
c_phi = 7500;                    %Nm/rad
d1 = 2*sqrt(c_phi*m_L);
c_theta = 1000;
d2 = 2*sqrt(c_theta*m_L);
c1 = 40 * 1000;

opts_flight = odeset('Events', @guard_flight);
opts_stance = odeset('Events', @guard_stance);
opts_apex = odeset('Events', @guard_apex);
step_width=.001;
x_full = [];
t_full = [];

t_apex = x0(5)/g;
tspan = 0:step_width:T;

[tout, xout, te, xe, ie] = ode45(@(t,x)mode2st(t,x,a,b,T,toffset), tspan, x0, opts_flight);
x_full = [x_full; xout];
t_full = [t_full; tout];


while isempty(xe) == false
    tspan = te:step_width:te+2*T;
    offset = xe(1) + L0 * sin(xe(3));
    x0 = xe';

    [tout, xout, te, xe, ie] = ode45(@(t,x)mode1st(t,x,a,b,T,offset,toffset), tspan, x0, opts_stance);
    x_full = [x_full; xout];
    t_full = [t_full; tout];
    
    
    if isempty(te) || sum(size(xe)~=[1,6])
        break
    end
    t_apex = te+xe(5)/g;
    tspan = te:step_width:te+2*T;
    x0 = xe';

    [tout, xout, te, xe, ie] = ode45(@(t,x)mode2st(t,x,a,b,T,toffset), tspan, x0, opts_flight);
    x_full = [x_full; xout];
    t_full = [t_full; tout];
    steps = steps +1;
end

figure; plot(x_full(:,1),x_full(:,2))
% figure; plot(x_full(:,2),x_full(:,5))
disp(strcat(int2str(steps),' steps made'));
