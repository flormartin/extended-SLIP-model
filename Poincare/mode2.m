function dx = mode2(t,x)
%x = [x y phi theta dx dy dphi dtheta]
% position x and y at center of entire Mass M along leg
% phi: upper leg angle / theta: lower leg angle

global L0 l0 l1 M m_B m_L g c_phi d1 c_theta d1 d2 c1 t_apex alpha0 omega;

%leg retraction phi0 = α0 + ω · (t − tapex)
phi_0 = alpha0 + omega * (t-t_apex);     %t_apex still missing

dx = [x(4);
    x(5);
    x(6);
    0.;
    -g;
    (-c_phi*(x(3)-phi_0)-d1*x(6)-m_L*g*sin(x(3))*l1)/(m_L*l1^2)];


end