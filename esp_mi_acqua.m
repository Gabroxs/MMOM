clearvars
close all
clc

short_struct=load('GammaShort.dat');
open_struct=load('GammaOpen.dat');

lshort=4.81e-3;
lopen=8e-3;
c = 3e8;
eps0 = 8.8541878176e-012;
mu0  = 1.2566370614e-006;

freq = short_struct(:,2) .* 1e9;
S11_open = db2mag(open_struct(:,3)) .* exp(1i .* deg2rad(open_struct(:,4)));
S11_short = db2mag(short_struct(:,3)) .* exp(1i .* deg2rad(short_struct(:,4)));
beta_aria = 2*pi .* freq .* sqrt(eps0 * mu0);

eps_r = ((1 - S11_open)./(1 + S11_open)) ./ (1i .* beta_aria * lopen);
mi_r = ((1 + S11_short)./(1 - S11_short)) ./ (1i .* beta_aria * lshort);

figure(1)
yyaxis left
plot(freq, real(eps_r), 'LineWidth',1.2);
ylabel('F/m')
yyaxis right
plot(freq, imag(eps_r), 'LineWidth',1.2);
ylabel('F/m')
xlabel('Hz')
title('\epsilon_r water')
legend('\Re[\epsilon_r]', '\Im[\epsilon_r]');

figure(2)
yyaxis left
plot(freq, real(mi_r), 'LineWidth',1.2);
ylabel('H/m')
yyaxis right
plot(freq, imag(mi_r), 'LineWidth',1.2);
ylabel('H/m')
xlabel('Hz')
title('\mu_r water')
legend('\Re[\mu_r]', '\Im[\mu_r]');
