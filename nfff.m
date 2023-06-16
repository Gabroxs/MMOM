%% NF FF transformation

clearvars
close all
clc

% Import data
load("campo2_MISURE_CORSO_MMM.mat");

%Plot campo misurato in z=d_misura
figure()
imagesc(x,y, 20*log10(abs(campo)/max(max(abs(campo)))))
xlabel('x') %misure in cm
ylabel('y') 
colorbar();
title('Campo misurato in z=d')

%% Parametri
c = 3e8;
omega = 2*pi*freq;
beta = (2*pi)/lambda;

%Passo di campionamento spaziale uniforme
dx = lambda/2;
dy = dx;

%Numero punti vettori x ed y
Nx = length(x);
Ny = length(y);

%Passo di campionamento nel dominio della frequenza
du = lambda/((Nx - 1)*dx);
dv = lambda/((Ny - 1)*dy);

u = -1:du:1;
v = -1:dv:1;

%Dominio del visibile

[U, V] = meshgrid(u,v);
domain = zeros(size(U));
domain(U.^2 + V.^2 <= 1) = 1;

%kz normalizzato rispetto k0
K = sqrt(1 - U.^2 - V.^2);

%Valuto lo spettro traslato in z=d
spec_d = fftshift(fft2(campo));

%Lo spettro in z=0 lo ottengo tramite lo shift
shift = domain .* exp(1i * beta * d_misura .* K);
spec_0 = spec_d .* shift;

figure()
imagesc(u,v, 20*log10(abs(spec_0)/max(max(abs(spec_0)))));
colorbar;
xlabel('u')
ylabel('v')
title('Spettro campo in z=0');


%Valuto il campo in z=0 antitrasformando lo spettro
campo_0 = ifft2(spec_0);

figure()
imagesc(x,y, 20*log10(abs(campo_0)/max(max(abs(campo_0)))));
xlabel('x')
ylabel('y')
title('Campo in z=0');
colorbar;


%Moltiplico il campo per la maschera
[X, Y] = meshgrid(x,y);
mask = zeros(size(X));
mask(25:38, 25:38) = 1;

figure()
campo_0_filtered = campo_0 .* mask;
imagesc(x,y, 20*log10(abs(campo_0_filtered)/max(max(abs(campo_0_filtered)))));
colorbar;
xlabel('x')
ylabel('y')
title('Campo in z=0 filtrato con maschera')

%Valuto lo spettro del campo filtrato in z=0 
spec_0_filtered = fft2(campo_0_filtered);
figure()
imagesc(u,v, 20*log10(abs(spec_0_filtered)/max(max(abs(spec_0_filtered)))));
colorbar;
xlabel('u')
ylabel('v')
title('Spettro campo filtrato in z=0')


%Far Field (spectrum) FF
ff = K.*spec_0_filtered;
figure()
imagesc(u,v, 20*log10(abs(ff)/max(max(abs(ff)))));
colorbar;
xlabel('u')
ylabel('v')
title('Far Field')


%Campo in FF
cff = ifft2(ff);
figure()
imagesc(x,y, 20*log10(abs(cff)/max(max(abs(cff)))));
colorbar;
xlabel('x')
ylabel('y')
title('Campo Far Field')