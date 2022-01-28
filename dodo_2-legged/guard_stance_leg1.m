function [value,isterminal,direction] = guard_stance_leg1(t,x)
%x = [x y phi theta dx dy dphi dtheta]
%Leg 1 is the leg with phi-angle

global offset l0 l1;
l = sqrt((offset-(x(1)+l1*sin(x(3))))^2+(x(2)-l1*cos(x(3)))^2);


value = l-l0;
isterminal = 1;     %stop the integration
direction = 1;      %from below

end
