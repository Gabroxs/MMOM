clearvars
close all
clc

%% Data import section

% Load data
temp = load('bb_mf.txt');
freq = temp(:,2);
bb_m = [temp(:,3), temp(:,4)];
bb_f = [temp(:,5), temp(:,6)];

% Open data
temp = load('open_mf.txt');
open_m = [temp(:,3), temp(:,4)];
open_f = [temp(:,5), temp(:,6)];

% Short data
temp = load('short_mf.txt');
short_m = [temp(:,3), temp(:,4)];
short_f = [temp(:,5), temp(:,6)];


% Through data
temp = load('through.txt');
through_11 = [temp(:,3), temp(:,4)];
through_12 = [temp(:,5), temp(:,6)];
through_21 = [temp(:,7), temp(:,8)];
through_22 = [temp(:,9), temp(:,10)];


% P1P3 no cal
temp = load('P1P3_no_cal.txt');
S11_no_cal = [temp(:,3) temp(:,4)];
S13_no_cal = [temp(:,5) temp(:,6)];
S31_no_cal = [temp(:,7) temp(:,8)];
S33_no_cal = [temp(:,9) temp(:,10)];
clear temp 

% P1P3 cal
temp = load('P1P3_cal.txt');
S11_cal = [temp(:,3) temp(:,4)];
S13_cal = [temp(:,5) temp(:,6)];
S31_cal = [temp(:,7) temp(:,8)];
S33_cal = [temp(:,9) temp(:,10)];
clear temp 

freq = freq * 1e9;
N = length(freq);
%% dB to complex conversion

% Load dB to complex
bb_f_cpx = db2mag(bb_f(:,1)).* exp(1i.*deg2rad(bb_f(:,2)));
bb_m_cpx = db2mag(bb_m(:,1)).* exp(1i.*deg2rad(bb_m(:,2)));

% Open dB to complex
open_f_cpx = db2mag(open_f(:,1)).* exp(1i.*deg2rad(open_f(:,2)));
open_m_cpx = db2mag(open_m(:,1)).* exp(1i.*deg2rad(open_m(:,2)));

% Short dB to complex
short_f_cpx = db2mag(short_f(:,1)).* exp(1i.*deg2rad(short_f(:,2)));
short_m_cpx = db2mag(short_m(:,1)).* exp(1i.*deg2rad(short_m(:,2)));

% Through dB to complex
through_11_cpx = db2mag(through_11(:,1)).* exp(1i.*deg2rad(through_11(:,2)));
through_12_cpx = db2mag(through_12(:,1)).* exp(1i.*deg2rad(through_12(:,2)));
through_21_cpx = db2mag(through_21(:,1)).* exp(1i.*deg2rad(through_21(:,2)));
through_22_cpx = db2mag(through_22(:,1)).* exp(1i.*deg2rad(through_22(:,2)));

%P1P3nocal dB to complex
S11_no_cal_cpx = db2mag(S11_no_cal(:,1)).* exp(1i.*deg2rad(S11_no_cal(:,2)));
S13_no_cal_cpx = db2mag(S13_no_cal(:,1)).* exp(1i.*deg2rad(S13_no_cal(:,2)));
S31_no_cal_cpx = db2mag(S31_no_cal(:,1)).* exp(1i.*deg2rad(S31_no_cal(:,2)));
S33_no_cal_cpx = db2mag(S33_no_cal(:,1)).* exp(1i.*deg2rad(S33_no_cal(:,2)));

%P1P3 cal dB to complex
S11_cal_cpx = db2mag(S11_cal(:,1)).* exp(1i.*deg2rad(S11_cal(:,2)));
S13_cal_cpx = db2mag(S13_cal(:,1)).* exp(1i.*deg2rad(S13_cal(:,2)));
S31_cal_cpx = db2mag(S31_cal(:,1)).* exp(1i.*deg2rad(S31_cal(:,2)));
S33_cal_cpx = db2mag(S33_cal(:,1)).* exp(1i.*deg2rad(S33_cal(:,2)));

%% Data and constants
Z0 = 50;
c = 3e8;
beta = (2*pi.*freq)/c;

% Gamma open
C0 = 16e-15;
C1 = -400e-27;
C2 = 35e-36;
C3 = 2.2e-45;
l_open = 8.966e-3;
C = C0 + C1*freq + C2*freq.^2 + C3*freq.^3;
Zc = 1./(1i*2*pi.*freq.*C);
Gamma_C = (Zc-Z0)./(Zc + Z0);
Gamma_open = Gamma_C .* exp(-2*beta*l_open*1i);
clear Gamma_C;

% Gamma short
l_short = 8.966e-3;
Gamma_S = -1;
Gamma_short = Gamma_S .* exp(-2*beta*l_short*1i);
clear Gamma_S;

% Gamma load
Gamma_load = zeros(N,1);

%% OSL calibration

% M matrix (forward)

M_osl_fw = zeros(3,3);
X_fw = zeros(3,1);

for j = 1:N

    M_osl_fw = [Gamma_open(j), 1, Gamma_open(j) * open_m_cpx(j);
                Gamma_short(j), 1, Gamma_short(j) * short_m_cpx(j);
                Gamma_load(j), 1, Gamma_load(j) * bb_m_cpx(j)    
                ];

    Gamma_m_fw = [open_m_cpx(j), short_m_cpx(j), bb_m_cpx(j)].';
    X_fw(:,j) = M_osl_fw\Gamma_m_fw;
        
end
clear j;

Er_f = X_fw(1,:);
Ed_f = X_fw(2,:);
Es_f = X_fw(3,:);


% M matrix (reverse)

M_osl_rev = zeros(3,3);
X_rev = zeros(3,N);

for j = 1:N

    M_osl_rev = [Gamma_open(j), 1, Gamma_open(j) * open_f_cpx(j);
                Gamma_short(j), 1, Gamma_short(j) * short_f_cpx(j);
                Gamma_load(j), 1, Gamma_load(j) * bb_f_cpx(j)    
                ];

    Gamma_m_rev = [open_f_cpx(j), short_f_cpx(j), bb_f_cpx(j)].';
    X_rev(:,j) = M_osl_rev\Gamma_m_rev;
    

end
clear j;

Er_r = X_rev(1,:);
Ed_r = X_rev(2,:);
Es_r = X_rev(3,:);


%% Through measurements

Ex_r = zeros(1,N);
Ex_f = zeros(1,N);

%forward and reverse
El_f = zeros(1,N);
Et_f = zeros(1,N);
El_r = zeros(1,N);
Et_r = zeros(1,N);

for j = 1:N
    %fw
    El_f(:,j) = (through_11_cpx(j,1) - Ed_f(1,j))/(through_11_cpx(j,1) * Es_f(1,j) + Er_f(1,j) - Es_f(1,j) * Ed_f(1,j));
    Et_f(:,j) = (through_21_cpx(j,1) - Ex_f(1,j)) * (1 - Es_f(1,j) * El_f(1,j));
    %rev
    El_r(:,j) = (through_22_cpx(j,1) - Ed_r(1,j))/(through_22_cpx(j,1) * Es_r(1,j) + Er_r(1,j) - Es_r(1,j) * Ed_r(1,j));
    Et_r(:,j) = (through_12_cpx(j,1) - Ex_r(1,j)) * (1 - Er_f(1,j) * El_r(1,j));
end

D = zeros(1,N);
S11 = zeros(1,N);
S13 = zeros(1,N);
S33 = zeros(1,N);
S31 = zeros(1,N);

for j = 1:N
    
    D(1,j) = (1 + ((S11_no_cal_cpx(j,1) - Ed_f(1,j))/Er_f(1,j))*Es_f(1,j)) * (1 + ((S33_no_cal_cpx(j,1) - Ed_r(1,j))*(Es_r(1,j)/Er_r(1,j)))) - ((S31_no_cal_cpx(j,1) - Ex_f(1,j))/Et_f(1,j))*((S31_no_cal_cpx(j,1) - Ex_r(1,j))/Et_r(1,j))*El_f(1,j)*El_r(1,j);
    S11(1,j) = (((S11_no_cal_cpx(j,1) - Ed_f(1,j))/Er_f(1,j)) * (1 + (Es_r(1,j)*(S33_no_cal_cpx(j,1) - Ed_r(1,j))/Et_r(1,j))) - El_f(1,j) * ((S13_no_cal_cpx(j,1) - Ex_r(1,j))/Et_r(1,j)) * ((S31_no_cal_cpx(j,1) - Ex_f(1,j))/Et_f(1,j)))/D(1,j);
    S13(1,j) = (((S31_no_cal_cpx(j,1) - Ex_f(1,j))/Et_f(1,j)) * (1 + ((Es_r(1,j)-El_f(1,j))*(S33_no_cal_cpx(j,1) - Ed_r(1,j))/Er_r(1,j))))/D(1,j);
    S33(1,j) = (((S33_no_cal_cpx(j,1) - Ed_r(1,j))/Er_r(1,j)) * (1 + Es_f(1,j) * (S11_no_cal_cpx(j,1) - Ed_f(1,j)/Er_f(1,j))) - El_r(1,j) * ((S31_no_cal_cpx(j,1) - Ex_f(1,j))/Et_f(1,j)) * ((S13_no_cal_cpx(j,1) - Ex_r(1,j))/Et_r(1,j)))/D(1,j); 
    S31(1,j) = (((S13_no_cal_cpx(j,1) - Ex_r(1,j))/Et_r(1,j)) * (1 + (Es_f(1,j) - El_r(1,j) )*((S11_no_cal_cpx(j,1) - Ed_f(1,j))/Er_f(1,j))))/D(1,j);
end 

figure()
plot(freq, 20*log10(abs(S11)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S11_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S11_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [GHz]')
ylabel('|S| [dB]')
legend('S11 matlab', 'S11 cal', 'S11 no cal')
title('S11 comparison')

figure()
plot(freq, 20*log10(abs(S13)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S13_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S13_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [GHz]')
ylabel('|S| [dB]')
legend('S13 matlab', 'S13 cal', 'S13 no cal')
title('S13 comparison')


figure()
plot(freq, 20*log10(abs(S33)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S33_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S33_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [GHz]')
ylabel('|S| [dB]')
legend('S33 matlab', 'S33 cal', 'S33 no cal')
title('S33 comparison')

figure()
plot(freq, 20*log10(abs(S31)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S31_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, 20*log10(abs(S31_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [GHz]')
ylabel('|S| [dB]')
legend('S31 matlab', 'S31 cal', 'S31 no cal')
title('S31 comparison')

figure()
plot(freq, rad2deg(angle(S31)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S31_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S31_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [Hz]')
ylabel('phase S [deg]')
legend('S31 matlab', 'S31 cal', 'S31 no cal')
title('S31 phase comparison')

figure()
plot(freq, rad2deg(angle(S11)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S11_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S11_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [Hz]')
ylabel('phase S [deg]')
legend('S11 matlab', 'S11 cal', 'S11 no cal')
title('S11 phase comparison')

figure()
plot(freq, rad2deg(angle(S33)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S33_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S33_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [Hz]')
ylabel('phase S [deg]')
legend('S33 matlab', 'S33 cal', 'S33 no cal')
title('S33 phase comparison')

figure()
plot(freq, rad2deg(angle(S13)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S13_cal_cpx)),'LineWidth', 1.15);
hold on
plot(freq, rad2deg(angle(S13_no_cal_cpx)), 'LineWidth', 0.8, 'LineStyle','--');
xlabel('frequency [Hz]')
ylabel('phase S [deg]')
legend('S13 matlab', 'S13 cal', 'S13 no cal')
title('S13 phase comparison')