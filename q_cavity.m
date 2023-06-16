% Misura Q-unloaded risonatore
clearvars
close all
clc

f0 = 4.9995e9;
f1 = 5.002e9;
delta_f = f1 - f0;
fbw = 2*delta_f/f0;
s21_ris = 708.613e-3;    
s21_f1 = 818.369e-3;    
g = (1 - s21_ris)/s21_ris;

Q = fbw^(-1) * sqrt((s21_f1^2 * (1 + g)^2 - 1)/(1 - s21_f1^2));
