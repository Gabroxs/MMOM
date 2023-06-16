clearvars 
close all
clc

%% Data import
temp = load("Guida_Banda_C_PLA_3_Spessore_5mm_Length_30.4_SParameter1.txt");
freq = temp(:,1);
S11 = temp(:,2) .* exp(1i*deg2rad(temp(:,3)));
S12 = temp(:,6) .* exp(1i*deg2rad(temp(:,7)));
S21 = temp(:,4) .* exp(1i*deg2rad(temp(:,5)));
S22 = temp(:,8) .* exp(1i*deg2rad(temp(:,9)));
%clear temp

%% Constants 
c = 3e8;
a = 34.8488e-3;
omega = 2*pi .* freq;
l = 30.4e-2;
l_mut = 5e-3;
l_shift = (l/2) - (l_mut/2);
lambda_0 = c./freq;
lambda_c = 2*a;
kz = 2*pi*sqrt((1./lambda_0).^2 - (1./lambda_c)^2);
theta = l_shift .* kz;
%% Gamma and T evaluation

%S_shift = diag(exp(-1i.*theta));

S11_sh = S11 .* exp(2i*theta);
S21_sh = S21 .* exp(2i*theta);

% S11_sh = exp(-1i.*theta) .* S11;
% S12_sh = exp(-1i.*theta) .* S12;
% S21_sh = exp(-1i.*theta) .* S21;
% S22_sh = exp(-1i.*theta) .* S22;

b = (S11_sh.^2 - S21_sh.^2 + 1)./(2.*S11_sh);

Gamma = b - sqrt(b.^2 - 1);

if (abs(Gamma) > 1)
    
    Gamma = b + sqrt(b.^2 - 1);

end

T = (S11_sh + S21_sh - Gamma)./(1 - (S11_sh + S21_sh).*Gamma);

% LL = (-1i/(2*pi*l_mut)) .* log(1./T);
% LL = (-1i/(2*pi*l_mut)) .* (log(1./(abs(T))) + 1i.*imag(1./T));
% 
% mi_r = LL .* ((1 + Gamma)./(1 - Gamma)) .* (1./sqrt((1./lambda_0).^2 - (1/lambda_c)^2));
% eps_r = ((1/lambda_c)^2 - (log(1./T)./(2*pi*l_mut)).^2) .*(lambda_0./mi_r).^2;
% %eps_r = ((1/lambda_c)^2 - (log(1./(abs(T))) + 1i.*imag(T))./2*pi*l_mut).^2 .*(lambda_0./mi_r).^2;

Log = log(1./(abs(T))) + 1i.*angle(1./T);
Lambda = (2 * pi * 1i * l_mut) ./ Log;
mi_r = ((1 + Gamma)./(Lambda .* (1 - Gamma))) ./ sqrt((1./lambda_0).^2 - (1./lambda_c).^2);
eps_r = (lambda_0 ./ mi_r).^2 .* ((1/lambda_c)^2 - (Log ./ (2 * pi * l_mut)).^2);

figure(1)
yyaxis left
plot(freq, real(eps_r),'LineWidth',1.2);
ylabel('F/m')
xlabel('Hz')
yyaxis right
plot(freq, imag(eps_r), 'LineWidth',1.2);
ylabel('F/m')
title('\epsilon_r')
legend('\Re[\epsilon_r]', '\Im[\epsilon_r]')

figure(2)
yyaxis left
plot(freq, real(mi_r),'LineWidth',1.2);
ylabel('H/m')
yyaxis right
plot(freq, imag(mi_r), 'LineWidth',1.2);
ylabel('H/m')
xlabel('Hz')
title('\mu_r ')
legend('\Re[\mu_r]', '\Im[\mu_r]')


% Gamma1 = b - sqrt(b.^2 - 1);
% Gamma2 = b + sqrt(b.^2 - 1);
% 
% abs(Gamma1)
% abs(Gamma2)


