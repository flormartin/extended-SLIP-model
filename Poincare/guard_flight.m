function [value,isterminal,direction] = guard_flight(t,x)
%
global L0;

value = x(2) - L0* cos(x(3));
isterminal = 1; %stop the integration
direction = -1; %from above

end

