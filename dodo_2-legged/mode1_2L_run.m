function dx = mode1_2L_run(t, x)
% x = [x y phi phi1 theta theta1 dx dy dphi dphi1 dtheta dtheta1]
% position x and y at center of entire Mass M along leg
% phi: angle of standing-foot / theta: angle of retracting foot

% import parameters
global g M m_L m_B l0 l1 c_hip c_knee dhip dknee c1 offset phi0 theta0 omega alpha0 t_apex;

%leg retraction phi0 = α0 + ω · (t − tapex)
u = alpha0 + omega * (t-t_apex);     %t_apex still missing

%phi0 = -35*pi/180;                              %motor set position
l = sqrt((offset-(x(1)+l1*sin(x(3))))^2+(x(2)-l1*cos(x(3)))^2);  %current spring length

%phi1 = x(3) - atan((offset-(x(1)+sin(x(3))*l1))/(x(2)-cos(x(3))*l1));
%theta1 = x(3) - atan((offset-(x(1)+sin(x(3))*l1))/(x(2)-cos(x(3))*l1));
%dtheta = x(6)-((-l1*cos(x(3))*x(6)-x(4))/(x(2)-l1*cos(x(3)))-((-l1*sin(x(3))+offset-x(1))*(l1*sin(x(3))*x(6)+x(5)))/(x(2)-l1*cos(x(3)))^2)/((-l1*sin(x(3))+offset-x(1))^2/(x(2)-l1*cos(x(3)))^2+1);


ddx =-c1*(l0-l)*sin(x(3)-x(4))/M;
ddy =-g + c1*(l0-l)*cos(x(3)-x(4))/M;
ddphi = (- m_L*g*l1*sin(x(3)) + m_L*g*l1*sin(x(5)) - x(9)*dhip - (x(4)-phi0)*c_hip + c1*(l0-l)*cos(x(9)-x(10))*l1*sin(x(9)) - c1*(l0-l)*sin(x(9)-x(10))*l1*cos(x(9)) - c_knee*x(4) - dknee*x(10) ) /(m_L*l1^2);
ddtheta = (- m_L*g*l1*sin(x(3)) + m_L*g*l1*sin(x(5)) - dhip*x(11) - c_hip*(x(5)-theta0) ) /(m_L*l1^2);

% derivative of state vector
dx = [x(7);
    x(8);
    x(9);
    x(9)-((-l1*cos(x(3))*x(9)-x(7))/(x(2)-l1*cos(x(3)))-((-l1*sin(x(3))+offset-x(1))*(l1*sin(x(3))*x(9)+x(8)))/(x(2)-l1*cos(x(3)))^2)/((-l1*sin(x(3))+offset-x(1))^2/(x(2)-l1*cos(x(3)))^2+1);
    x(11);
    0; %%no angular speed on spring-to-knee for retracting foot
    ddx;
    ddy;
    ddphi;
    (((2*(-l1*sin(x(3)) + offset - x(1))*(-l1*cos(x(3))*x(9) - x(7)))/(x(2) - l1*cos(x(3)))^2 - (2*(-l1*sin(x(3)) + offset - x(1))^2 *(l1*sin(x(3))*x(9) + x(8)))/(x(2) - l1*cos(x(3)))^3)*((-l1*cos(x(3))*x(9) - x(7))/(x(2) - l1*cos(x(3))) - ((-l1*sin(x(3)) + offset - x(1))*(l1*sin(x(3))*x(9) + x(8)))/(x(2) - l1*cos(x(3)))^2))/((-l1*sin(x(3)) + offset - x(1))^2/(x(2) - l1*cos(x(3)))^2 + 1)^2 - (-((-l1*sin(x(3)) + offset - x(1))*(l1*cos(x(3))*x(9)^2 + l1*sin(x(3))*ddphi + ddy))/(x(2) - l1*cos(x(3)))^2 + (2*(-l1*sin(x(3)) + offset - x(1))*(l1*sin(x(3))*x(9) + x(8))^2)/(x(2) - l1*cos(x(3)))^3 + (l1*sin(x(3))*x(9)^2 - l1*cos(x(3))*ddphi - ddx)/(x(2) - l1*cos(x(3))) - (2*(-l1*cos(x(3))*x(9) - x(7))*(l1*sin(x(3))*x(9) + x(8)))/(x(2) - l1*cos(x(3)))^2)/((-l1*sin(x(3)) + offset - x(1))^2/(x(2) - l1*cos(x(3)))^2 + 1) + ddphi;
    ddtheta + u; %%retracting foot is actuated
    0]; %%no angular acceleration on spring-to-knee for retracting foot
end
