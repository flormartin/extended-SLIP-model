function dx = mode2(t,x)
%x = [x y phi theta dx dy dphi dtheta]
% position x and y at center of entire Mass M along leg
% phi: upper leg angle / theta: lower leg angle

M = 80.;
m_leg = .32*M;
m_trunk = M-m_leg;
L0 = 1;
k_hip = ; %set value
d_hip = ; %set value
omega_ret = 50*pi/180;
alpha0 = 68 * pi / 180;

%leg retraction phi0 = α0 + ω · (t − tapex)
phi0 = alpha0 + omega_ret * (t-t_apex);     %t_apex still missing

dx = [x(5);                                                 %dx
    x(6);                                                   %dy
    x(7);                                                   %dphi
    x(8);                                                   %dtheta
    0.;                                                     %ddx
    -9.81;                                                  %ddy
    (k_hip*(x(3)-phi0)+d_hip*(omega_ret-x(7)))/...          %ddphi
    (m_leg*(m_trunk/M*.4*L0)^2+m_trunk*(mleg/M*.4*L0)^2);                    %
    0.];                                                    %ddtheta
end