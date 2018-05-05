%% Cross-Frequency Coupling

%% Hilbert transform

% Generate sine wave
dt = 0.001; % time step [s]
t = dt:dt:1; % time axis [s]
f1 = 2.0; % freq of sine wave [Hz]
phi0 = 0.0; % initial phase of sinusoid
d = sin(2.0*pi*t*f1 + phi0); 

% Compute analytic signal (taking the Hilbert transform) 
dA = hilbert(d); % complex = real + imag

% Plot
figure;
subplot(4,1,1); plot(d); ylabel('Data');
subplot(4,1,2); plot(real(dA)); 
hold on;
plot(imag(dA), 'r'); 
hold off;
ylabel('Real (blue), Imag (red)');
axis tight

% Compute phase
phi = angle(dA);

% Compute amplitude envelope
amp = abs(dA);

% Plot results 
subplot(4,1,3); plot(phi); ylabel('Angle'); axis tight
subplot(4,1,4); plot(amp); ylabel('Amplitude'); axis tight

%% Power Spectrum
% units: [watts], [V^2]
% power density: [watts]/[Hz] or [V^2]/[Hz]

% Generate a noisy signal
dt = 0.001; 
t = dt:dt:1; 
T = max(t);
N = length(t);
d = randn(N,1); % normalized 

% Compute analytic signal
dA = hilbert(d);

pow_d = 2*dt^2/T * fft(d).*conj(fft(d)); 
pow_dA = 2*dt^2/T * fft(dA).*conj(fft(dA)); 

% Define frequency axis
df = 1/T;

% Nyquist frequency: If you have a sampling rate of "v," "fNQ" = 0.5v

fNQ = 1/dt/2;
faxis = -fNQ: df: fNQ-df;

figure;
subplot(2,1,1); plot(t,d); ylabel('Data'); xlabel('Time [s]');
subplot(2,1,2); plot(faxis, pow_d); ylabel('Power of Signal'); xlabel('Freq [Hz]');
hold on;
plot(faxis, pow_dA/4, 'r'); title('Power of Signal (red)'); hold off; 

%% Compute CFC
load('data_1.mat'); 

% plot data and check for CFC by eye
t = (1:length(d))*dt; % time axis
N = length(t);
Fs = 1/dt; % sampling freq
plot(t, d); 
xlabel('Time [s]');

periodogram(d, [], [], 1/dt); 

% Filter the data into the low and high freq bands identified by eye
% Butterworth filter, second order

deg = 2; % filter order
Wn = [5*2/Fs, 7*2/Fs]; % low freq window of interest
[B,A] = butter(deg, Wn, 'bandpass'); % apply filter to isolate band
dlo = filtfilt(B,A,d); 

Wn = [60*2/Fs, 140*2/Fs]; % high freq window of interest
[B,A] = butter(deg, Wn, 'bandpass');
dhi = filtfilt(B,A,d);

% How well does the filter work?
plot(t,d);
hold on;
plot(t, dlo, 'r');
plot(t, dhi, 'g');
hold off;
xlabel('Time [s]'); 
title('Data (blue), Low-freq (red), High-freq (green)');














































































































