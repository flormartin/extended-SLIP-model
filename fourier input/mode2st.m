function dx = mode2(t,x,a,b,T)
%x = [x y phi theta dx dy dphi dtheta]
% position x and y at center of entire Mass M along leg
% phi: upper leg angle / theta: lower leg angle

phi0 = a(1) + sum(a(2).*cos(t*2*pi/T)+b.*sin(t*2*pi/T),1);

global g m_L l1 c_phi d1;

dx = [x(4);
    x(5);
    x(6);
    0.;
    -g;
    (-c_phi*(x(3)-phi0)-d1*x(6)-m_L*g*sin(x(3))*l1)/(m_L*l1^2)];
end