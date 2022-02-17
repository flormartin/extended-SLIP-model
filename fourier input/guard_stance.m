function [value,isterminal,direction] = guard_stance(t,x)
%x = [x y phi dx dy dphi]

global offset l0 l1;
x_knee = x(1) + l1 * sin(x(3)) - offset; 
y_knee = x(2) - l1 * cos(x(3));
l = sqrt(x_knee^2+y_knee^2);

value = l-l0;       %liftoff
isterminal = 1;     %stop the integration
direction = 1;      %from below

end