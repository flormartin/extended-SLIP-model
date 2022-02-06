clear;


global L0 l0 l1 M m_B m_L g c_phi d1 c_theta d2 c1 t_apex alpha0 omega;
L0=1; M=80; g=9.81;
l0 = .6*L0;                     %lower leg length
l1 = .4*L0;                     %upper leg length
m_L = .32*M;                    %leg mass
m_B = M-m_L;                    %rest of body mass
c_phi = 750;                    %Nm/rad
d1 = 2*sqrt(c_phi*m_L);
c_theta = 1000;
d2 = 2*sqrt(c_theta*m_L);
c1 = 20*1000;
%alpha0 = (90 - 62) * pi/ 180.;
% t_apex = X(5,1)/g;              %vertical initial speed divided by g
%omega = 50*pi/180; %rad/s

alpha0 = (90 - 62) * pi/ 180.;
t_apex = 0;
omega = 50*pi/180; %rad/s


% initial guess
X0 = [0.; 1.05; -35 * pi/ 180.; 2.5; 2; 0];
fprintf('eig(monodromy(X0)) = %f\n',eig(monodromy(X0)));

% X0 = [0.; .8; -35 * pi/ 180.; 5; .5; 0];
y_guess = X0(2);

X_old = X0;
y_old = X0(2);

error = 1000.;
while error > 0.001
    
[X_pred, X_full, te] = poincare(X_old);
y_pred = X_pred(2);
y_update = (1./(1. - monodromy(X_old))) * (y_pred - y_old);

y_next = y_old + y_update;

X_old(2) = y_next;
y_old = y_next;

error = abs(y_old-y_pred);

end

plot(X_full(:,1), X_full(:,2))





