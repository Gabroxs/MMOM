clc; close all; clear;
% file_dati = 'S11CavoNero.txt';
% [f, S11] = S11Extract(file_dati);

temp = load("S11CavoNeroTIMEHARMONIC_mod.txt");
freq = temp(:,2) .* 1e9;
N = length(freq);
% S11_db = temp(:,3);
% S11_deg = temp(:,4);
S11 = db2mag(temp(:,3)) .* exp(1i .* deg2rad(temp(:,4)));

f_ext = [-flip(freq);0;freq];
S11_ext = [conj(flip(S11));S11(1);S11];

%plot(f,(20*log10(S11)));
%figure
%plot(f_ext,(20*log10(S11_ext)));

%Parameters
c = 3e8;
vp = 0.82 * c;
%pulse_width = 1/(freq(end) - freq(1));
%delta_f = (freq(end) - freq(1))/(N - 1);
delta_f = freq(end) - freq(1);
max_time_span = 1/(2*delta_f);

band = 4e9;
t = (0:N-1)*max_time_span;
S11_TD = ifft(S11_ext);
RL_TD = 20*log10(abs(S11_TD));
l = vp*t/2;

figure()
plot(l,RL_TD(1:N),'LineWidth',1.15);
xlabel('distance [m]')
ylabel('|S11| [dB]')
[max,index] = max(RL_TD);
l_cavo = l(index)