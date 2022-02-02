function dx = stance(x,phi0,offset)
% x = [x y phi theta dx dy dphi dtheta]
% position x and y at center of entire Mass M along leg
% phi: upper leg angle / theta: lower leg angle

% import parameters
global g M m_L l0 l1 c_phi c_theta d1 d2 c1;

%phi0 = -35*pi/180;                              %motor set position
l = sqrt((offset-(x(1)+l1*sin(x(3))))^2+(x(2)-l1*cos(x(3)))^2);  %current spring length

theta = x(3) - atan((offset-(x(1)+sin(x(3))*l1))/(x(2)-cos(x(3))*l1));
dtheta = x(6)-((-l1*cos(x(3))*x(6)-x(4))/(x(2)-l1*cos(x(3)))-((-l1*sin(x(3))+offset-x(1))*(l1*sin(x(3))*x(6)+x(5)))/(x(2)-l1*cos(x(3)))^2)/((-l1*sin(x(3))+offset-x(1))^2/(x(2)-l1*cos(x(3)))^2+1);

ddx =-c1*(l0-l)*sin(x(3)-theta)/M;
ddy =-g + c1*(l0-l)*cos(x(3)-theta)/M;
ddphi = (-c_phi*(x(3)-phi0)     -d1*x(6) -c_theta*theta   -d2*dtheta    -m_L*g*sin(x(3))*l1+c1*(l0-l)*sin(x(3)-theta)*l1)    /(m_L*l1^2);

% derivative of state vector
dx = [x(4);
    x(5);
    x(6);
    ddx;
    ddy;
    ddphi];
end