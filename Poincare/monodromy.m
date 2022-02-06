function dP = monodromy(x)

r = zeros(6,1);
epsilon = 1e-5;
r(2) = epsilon;

df =  (poincare(x+r) - poincare(x-r)) / (2.*epsilon);
dP = df(2);

end

