%%Load 
z1 = load('\\zserver.cortexlab.net\Code\Rigging\ExpDefinitions\Pip\valveTestZYM1.mat');
z2 = load('\\zserver.cortexlab.net\Code\Rigging\ExpDefinitions\Pip\valveTestZYM2.mat');

Z1 = cat(1,z1.audD{:});
Z2 = cat(1,z2.audD{:});

%%

Fs = 192e3;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(Z1);             % Length of signal
t = (0:L-1)*T;        % Time vector


figure;
subplot(2,2,1);
plot(t,Z1); title('Zym1');
subplot(2,2,2);
plot(t,Z2); title('Zym2');



Y = fft(Z1);
P2 = abs(Y/L);
P1a = P2(1:L/2+1);
P1a(2:end-1) = 2*P1a(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,3);
plot(f,P1a);  
title('Spectrum for Zym1')
xlabel('f (Hz)')
ylabel('|P1(f)|')

Y = fft(Z2);
P2 = abs(Y/L);
P1b = P2(1:L/2+1);
P1b(2:end-1) = 2*P1b(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,4);
plot(f,P1b);  
title('Spectrum for Zym1')
xlabel('f (Hz)')
ylabel('|P1(f)|')