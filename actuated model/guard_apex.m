function [value,isterminal,direction] = guard_apex(t,x)
%x = [x y phi dx dy dphi]

value = x(5);   %apex during flight phase
isterminal = 1; %stop the integration
direction = -1; %from above

end

