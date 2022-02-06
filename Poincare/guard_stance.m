function [value,isterminal,direction] = guard_stance(t,x)
%x = [x y dx dy]

% global offset l0 l1;
global L0 l0 l1;

phi = x(3);
offset = x(1) + L0 * sin(phi);
l = sqrt((offset-(x(1)+l1*sin(x(3))))^2+(x(2)-l1*cos(x(3)))^2);


value = l-l0;
isterminal = 1;     %stop the integration
direction = 1;      %from below

end