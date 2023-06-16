clearvars 
close all
clc

%% Data import section

%load = 1;  %load = 1 -> gamma_l = -1; load = 0 -> gamma_l = 1

%through data
S11_T1 = load('thruS11nocal.txt');
S12_T1 = load('thruS12nocal.txt');
S21_T1 = load('thruS21nocal.txt');
S22_T1 = load('thruS22nocal.txt');

%line data
S11_T2 = load('lineS11nocal.txt');
S12_T2 = load('lineS12nocal.txt');
S21_T2 = load('lineS21nocal.txt');
S22_T2 = load('lineS22nocal.txt');

%Rm (reflect data)
S11_Rm = load("RXMS11nocal.txt");
S22_Rm = load('RYMS22nocal.txt');

freq = S11_T1(:,2) .* 1e9; %[Hz]

%isolator data (nocal)
S11_Tm_nocal = load('dutisolatoreS11nocal.txt');
S12_Tm_nocal = load('dutisolatoreS12nocal.txt');
S21_Tm_nocal = load('dutisolatoreS21nocal.txt');
S22_Tm_nocal = load('dutisolatoreS22nocal.txt');

%isolator data (cal)
S11_Tm_cal = load('dutisolatoreS11cal.txt');
S12_Tm_cal = load('dutisolatoreS12cal.txt');
S21_Tm_cal = load('dutisolatoreS21cal.txt');
S22_Tm_cal = load('dutisolatoreS22cal.txt');

%% Data conversion db,deg -> natural,rad

%through data
S11_T1 = db2mag(S11_T1(:,3)) .* exp(1i .* deg2rad(S11_T1(:,4)));
S12_T1 = db2mag(S12_T1(:,3)) .* exp(1i .* deg2rad(S12_T1(:,4)));
S21_T1 = db2mag(S21_T1(:,3)) .* exp(1i .* deg2rad(S21_T1(:,4)));
S22_T1 = db2mag(S22_T1(:,3)) .* exp(1i .* deg2rad(S22_T1(:,4)));

%line data
S11_T2 = db2mag(S11_T2(:,3)) .* exp(1i .* deg2rad(S11_T2(:,4)));
S12_T2 = db2mag(S12_T2(:,3)) .* exp(1i .* deg2rad(S12_T2(:,4)));
S21_T2 = db2mag(S21_T2(:,3)) .* exp(1i .* deg2rad(S21_T2(:,4)));
S22_T2 = db2mag(S22_T2(:,3)) .* exp(1i .* deg2rad(S22_T2(:,4)));

%Reflect
S11_Rm = db2mag(S11_Rm(:,3)) .* exp(1i .* deg2rad(S11_Rm(:,4)));
S22_Rm = db2mag(S22_Rm(:,3)) .* exp(1i .* deg2rad(S22_Rm(:,4)));

%isolator data (nocal)
S11_Tm_nocal = db2mag(S11_Tm_nocal(:,3)) .* exp(1i .* deg2rad(S11_Tm_nocal(:,4)));
S12_Tm_nocal = db2mag(S12_Tm_nocal(:,3)) .* exp(1i .* deg2rad(S12_Tm_nocal(:,4)));
S21_Tm_nocal = db2mag(S21_Tm_nocal(:,3)) .* exp(1i .* deg2rad(S21_Tm_nocal(:,4)));
S22_Tm_nocal = db2mag(S22_Tm_nocal(:,3)) .* exp(1i .* deg2rad(S22_Tm_nocal(:,4)));

%isolator data (cal)
S11_Tm_cal = db2mag(S11_Tm_cal(:,3)) .* exp(1i .* deg2rad(S11_Tm_cal(:,4)));
S12_Tm_cal = db2mag(S12_Tm_cal(:,3)) .* exp(1i .* deg2rad(S12_Tm_cal(:,4)));
S21_Tm_cal = db2mag(S21_Tm_cal(:,3)) .* exp(1i .* deg2rad(S21_Tm_cal(:,4)));
S22_Tm_cal = db2mag(S22_Tm_cal(:,3)) .* exp(1i .* deg2rad(S22_Tm_cal(:,4)));

%% TRL calibration
% T1 = S2T_fixed(S11_T1, S12_T1, S21_T1, S22_T1);     %through
% T2 = S2T_fixed(S11_T2, S12_T2, S21_T2, S22_T2);     %line
% T_nocal = S2T_fixed(S11_Tm_nocal, S12_Tm_nocal, S21_Tm_nocal, S22_Tm_nocal);    %dut no cal
% T_cal = S2T_fixed(S11_Tm_cal, S12_Tm_cal, S21_Tm_cal, S22_Tm_cal);              %dut cal

N = length(S11_T1);

for i = 1 : N

    T1 = 1/(S21_T1(i)) * [-((S11_T1(i) * S22_T1(i) - S12_T1(i) * S21_T1(i))), S11_T1(i);
                            -S22_T1(i), 1];

    T2 = 1/(S21_T2(i)) * [-((S11_T2(i) * S22_T2(i) - S12_T2(i) * S21_T2(i))), S11_T2(i);
                            -S22_T2(i), 1];

    T_dut_nocal = 1/(S21_Tm_nocal(i)) * [-((S11_Tm_nocal(i) * S22_Tm_nocal(i) - S12_Tm_nocal(i) * S21_Tm_nocal(i))), S11_Tm_nocal(i);
                            -S22_Tm_nocal(i), 1];

    T_dut_cal = 1/(S21_Tm_cal(i)) * [-((S11_Tm_cal(i) * S22_Tm_cal(i) - S12_Tm_cal(i) * S21_Tm_cal(i))), S11_Tm_cal(i);
                            -S22_Tm_cal(i), 1];


    T3 = T2/T1;
    delta = (T3(2,2) - T3(1,1))^2 + 4 * T3(1,2) * T3(2,1);

    sol_1 = (T3(1,1) - T3(2,2) - sqrt(delta))/(2 * T3(2,1));
    sol_2 = (T3(1,1) - T3(2,2) + sqrt(delta))/(2 * T3(2,1));

    if abs(sol_1) < abs(sol_2)

        bx = sol_1;
        ax_over_cx = sol_2;
    
    else

        bx = sol_2;
        ax_over_cx = sol_1;


    end




T1m = [T1(1,1)/T1(2,2) T1(1,2)/T1(2,2);
       T1(2,1)/T1(2,2) 1];

%alpha evaluation
alpha = (1 - T1m(1,2)/ax_over_cx)/(1 - bx/ax_over_cx);

% x22 * y22 evaluation
x22_y22 = alpha * T1(2,2);

% cy evaluation
cy = (T1m(2,1) - T1m(1,1)/ax_over_cx)/(alpha * (1 - bx/ax_over_cx));

% ax * ay evaluation
ax_ay = (T1m(1,1) - bx * T1m(2,1))/(alpha * (1 - bx/ax_over_cx));

%by/ay evaluation
by_over_ay = (T1m(1,2) - bx)/(T1m(1,1) - bx * T1m(2,1));

%ay/ax evaluation
ay_over_ax = ((cy + S22_Rm(i))/(1 + S22_Rm(i) * by_over_ay)) * ((1 - S11_Rm(i)/ax_over_cx) / (S11_Rm(i) - bx));

% A1 and A2 evaluation

A1 = ax_ay;
A2 = ay_over_ax;

% if load == 1
% 
%     gamma_l = -1;
% 
% else
% 
%     gamma_l = 1;
% 
% end
gamma_l = -1;

ay2 = (cy + S22_Rm(i))/(gamma_l * (1 + S22_Rm(i) * by_over_ay));

%ax and ay evaluation

sol_1 = sqrt(A1 * A2);
sol_2 = -sol_1;

phi_1 = angle(sol_1);
phi_2 = angle(sol_2);
phi = angle(ay2);

delta_phi_1 = phi - phi_1;
delta_phi_2 = phi - phi_2;

if abs(delta_phi_1) < abs(delta_phi_2)

    ay = sol_1;

else

    ay = sol_2;

end

ax = A1/ay;

%by and cy evaluation

by = by_over_ay * ay;
cx = ax/ax_over_cx;

%% T_dut matrix evaluation

T_dut = (inv([ax bx; cx 1]) * T_dut_nocal * inv([ay by; cy 1]))/x22_y22;
%S_dut = t2s(T_dut);
S_dut = t2s(rot90(T_dut,2));

S11_dut(i) = S_dut(1,1);
S12_dut(i) = S_dut(1,2);
S21_dut(i) = S_dut(2,1);
S22_dut(i) = S_dut(2,2);

end

figure()
plot(freq, 20*log10(abs(S11_dut)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S11_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S11_Tm_nocal)), 'LineStyle','--','LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('|S| [dB]');
legend('S11 Matlab', 'S11 cal', 'S11 no cal');
title('S11 isolator comparison')

figure()
plot(freq, 20*log10(abs(S12_dut)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S12_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S12_Tm_nocal)), 'LineStyle','--','LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('|S| [dB]');
legend('S12 Matlab', 'S12 cal', 'S12 no cal');
title('S12 isolator comparison')

figure()
plot(freq, 20*log10(abs(S21_dut)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S21_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S21_Tm_nocal)), 'LineStyle','--', 'LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('|S| [dB]');
legend('S21 Matlab', 'S21 cal', 'S21 no cal');
title('S21 isolator comparison')

figure()
plot(freq, 20*log10(abs(S22_dut)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S22_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, 20*log10(abs(S22_Tm_nocal)), 'LineStyle','--','LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('|S| [dB]');
legend('S22 Matlab', 'S22 cal', 'S22 no cal');
title('S22 isolator comparison')

figure()
plot(freq, rad2deg(angle(S11_dut)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S11_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S11_Tm_nocal)), 'LineStyle','--','LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('S phase [deg]');
legend('S11 Matlab', 'S11 cal', 'S11 no cal');
title('S11 isolator phase comparison')

figure()
plot(freq, rad2deg(angle(S12_dut)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S12_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S12_Tm_nocal)), 'LineStyle','--','LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('S phase [deg]');
legend('S12 Matlab', 'S12 cal', 'S12 no cal');
title('S12 isolator phase comparison')

figure()
plot(freq, rad2deg(angle(S21_dut)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S21_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S21_Tm_nocal)), 'LineStyle','--','LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('S phase [deg]');
legend('S21 Matlab', 'S21 cal', 'S21 no cal');
title('S21 isolator phase comparison')

figure()
plot(freq, rad2deg(angle(S22_dut)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S22_Tm_cal)), 'LineWidth',1.15);
hold on
plot(freq, rad2deg(angle(S22_Tm_nocal)), 'LineStyle','--','LineWidth',0.8);
xlabel('Frequency [Hz]');
ylabel('S phase [deg]');
legend('S22 Matlab', 'S22 cal', 'S22 no cal');
title('S22 isolator phase comparison')

