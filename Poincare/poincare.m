function [P, x_full, te] = poincare(x_previous)

global L0 l0 l1 M m_B m_L g c_phi d1 c_theta d2 c1 t_apex alpha0 omega;



x0 = x_previous;

opts_flight = odeset('Events', @guard_flight);
opts_stance = odeset('Events', @guard_stance);

step_width = 0.001;


x_full = [];
t_full = [];
%% 
% mode2 flight phase
% mode1 stance phase

tspan = 0:step_width:2.0;
[tout, xout, te, xe, ie] = ode45(@mode2, tspan, x0, opts_flight);
x_full = [x_full; xout];
t_full = [t_full; tout];


tspan = te:step_width:2.0;
phi = xe(3);
offset = xe(1) + L0 * sin(phi);

x0 = xe;
[tout, xout, te, xe, ie] = ode45(@mode1, tspan, x0,opts_stance);
x_full = [x_full; xout];
t_full = [t_full; tout];

% %% 
% tspan = te:step_width:2.0;
% [tout, xout, te, xe, ie] = ode45(@mode2, tspan, x0, opts_flight);
% x_full = [x_full; xout];
% t_full = [t_full; tout];
% 
% tspan = te:step_width:2.0;
% phi = xe(3);
% offset = xe(1) + L0 * sin(phi);
% 
% x0 = xe;
% [tout, xout, te, xe, ie] = ode45(@mode1, tspan, x0,opts_stance);
% x_full = [x_full; xout];
% t_full = [t_full; tout];

P = xe;

end
