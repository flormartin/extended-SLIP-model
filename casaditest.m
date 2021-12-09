import casadi.*

opti = casadi.Opti();

x = opti.variable();
y = opti.variable();

opti.minimize((1-x)^2+(y-x^2)^2);

opti.solver('ipopt');
sol = opti.solve();

plot(sol.value(x),sol.value(y),'o');