function dx = mode2_2L_run(t, x)
%x = [x y phi theta dx dy dphi dtheta]
% position x and y at center of entire Mass M along leg
% phi: angle of retracting-foot / theta: angle of standing foot

global g M m_L m_B l0 l1 c_hip c_knee dhip dknee theta0 c1 offset phi0 theta0 omega alpha0 t_apex;

%leg retraction phi0 = α0 + ω · (t − tapex)
u = alpha0 + omega * (t-t_apex);     %t_apex still missing

%theta0 = -35*pi/180;                              %motor set position
l = sqrt((offset-(x(1)+l1*sin(x(5))))^2+(x(2)-l1*cos(x(5)))^2);  %current spring length

ddx =-c1*(l0-l)*sin(x(5)-x(6))/M;
ddy =-g + c1*(l0-l)*cos(x(5)-x(6))/M;
ddphi = (- m_L*g*l1*sin(x(3)) + m_L*g*l1*sin(x(5)) - dhip*x(9) - c_hip*(x(3)-phi0) ) /(m_L*l1^2);
ddtheta = (- m_L*g*l1*sin(x(3)) + m_L*g*l1*sin(x(5)) - x(11)*dhip - (x(5)-theta0)*c_hip + c1*(l0-l)*cos(x(11)-x(12))*l1*sin(x(11)) - c1*(l0-l)*sin(x(11)-x(12))*l1*cos(x(11)) - c_knee*x(6) - dknee*x(12) ) /(m_L*l1^2);

dx = [x(7);
    x(8);
    x(9);
    0; %%no angular speed on spring-to-knee for retracting foot
    x(11);
    x(11)-((-l1*cos(x(5))*x(11)-x(7))/(x(2)-l1*cos(x(5)))-((-l1*sin(x(5))+offset-x(1))*(l1*sin(x(5))*x(11)+x(8)))/(x(2)-l1*cos(x(5)))^2)/((-l1*sin(x(5))+offset-x(1))^2/(x(2)-l1*cos(x(5)))^2+1);
    ddx;
    ddy;
    ddphi + u; %%retracting foot is actuated
    0; %%no angular acceleration on spring-to-knee for retracting foot
    ddtheta;
    (((2*(-l1*sin(x(5)) + offset - x(1))*(-l1*cos(x(5))*x(11) - x(7)))/(x(2) - l1*cos(x(5)))^2 - (2*(-l1*sin(x(5)) + offset - x(1))^2 *(l1*sin(x(5))*x(11) + x(8)))/(x(2) - l1*cos(x(5)))^3)*((-l1*cos(x(5))*x(11) - x(7))/(x(2) - l1*cos(x(5))) - ((-l1*sin(x(5)) + offset - x(1))*(l1*sin(x(5))*x(11) + x(8)))/(x(2) - l1*cos(x(5)))^2))/((-l1*sin(x(5)) + offset - x(1))^2/(x(2) - l1*cos(x(5)))^2 + 1)^2 - (-((-l1*sin(x(5)) + offset - x(1))*(l1*cos(x(5))*x(11)^2 + l1*sin(x(5))*ddtheta + ddy))/(x(2) - l1*cos(x(5)))^2 + (2*(-l1*sin(x(5)) + offset - x(1))*(l1*sin(x(5))*x(11) + x(8))^2)/(x(2) - l1*cos(x(5)))^3 + (l1*sin(x(5))*x(11)^2 - l1*cos(x(5))*ddtheta - ddx)/(x(2) - l1*cos(x(5))) - (2*(-l1*cos(x(5))*x(11) - x(7))*(l1*sin(x(5))*x(11) + x(8)))/(x(2) - l1*cos(x(5)))^2)/((-l1*sin(x(5)) + offset - x(1))^2/(x(2) - l1*cos(x(5)))^2 + 1) + ddtheta];
end