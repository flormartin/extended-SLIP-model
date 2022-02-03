clear all

N = 80;
n_coeff = 8;

N_phase = N/2;
a = randn(1,n_coeff+1);
b = randn(1,n_coeff);
k = 1:n_coeff;
T_fl=.25;
T_st=.15;
h_fl=T_fl/N_phase;
h_st=T_st/N_phase;
i = 25;
t = h_st * i;
phi0 = a(1) + sum(a(k+1).*cos(t*k/(T_st+T_fl))+b(k).*sin(t*k/(T_st+T_fl)))