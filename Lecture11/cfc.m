% Cross-Frequency Coupling Demo
% code written by Dr Mark Kramer, Boston University

%% Hilbert Transform

% Generate sinusoidal data
dt = 0.001; % sampling interval [s]
t = dt:dt:1; % time axis [s]
f1 = 2.0; % freq of sinusoid [Hz]
phi0 = 0.0; % initial phase of sinusoid
d = sin(2.0*pi*t*f1 + phi0);

% Compute analytic signal (complex-valued function with neg freq)
dA = hilbert(d); 

figure; 
subplot(4,1,1); plot(d); ylabel('Data');
subplot(4,1,2); plot(real(dA)); 
hold on;
plot(imag(dA),'r'); 
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
% Analytic signal has no power at negative frequencies. We will show this
% result in a numerical simulation

% Generate a noisy signal
dt = 0.001;
t = dt:dt:1;
T = max(t);
N = length(t); 
d = randn(N,1); 

% compute analytic signal
dA = hilbert(d); 

pow_d = 2*dt^2/T * fft(d).*conj(fft(d)); 
pow_dA = 2*dt^2/T * fft(dA).*conj(fft(dA)); 

% Define freqency axis
df = 1/T; 

% Nyquist frequency: In order to recover all Fourier components of a
% periodic waveform, it is necessary to use a sampling rate "v" at least
% twice the highest waveform frequency. so "fNQ" = 0.5v

fNQ = 1/dt/2;
faxis = -fNQ: df: fNQ-df; 

figure; 
subplot(2,1,1); plot(t,d); ylabel('Data'); xlabel('Time [s]');
subplot(2,1,2); plot(faxis, pow_d); ylabel('Power of Signal'); xlabel('Freq [Hz]'); 
hold on; 
plot(faxis, pow_dA/4, 'r'); title('Power of Signal (red)'); hold off; 

%% Compute CFC of signal
load('data_1.mat');

% plot data and check for CFC by eye
t = (1:length(d))*dt; % time axis
N = length(t); 
Fs = 1/dt; % sampling freq
plot(t, d); 
xlabel('Time [s]'); 

periodogram(d, [], [], 1/dt); 

% From this figure, there are two peaks:
% (i) sharp low freq peak near 6 Hz
% (ii) broad high freq peak from 60-140 Hz
% Note: high freq activity has much lower power than low freq activity

% Goal: develop measure to quantify CFC 
% Step 1: Filter data into low and high freq bands 

% we will use second order Butterworth filter

deg = 2; % filter order
Wn = [5*2/Fs, 7*2/Fs]; % low freq window of interest
[B,A] = butter(deg, Wn, 'bandpass'); % apply filter to isolate band
dlo = filtfilt(B,A,d); 

Wn = [60*2/Fs, 140*2/Fs]; % high freq window of interest
[B,A] = butter(deg, Wn, 'bandpass'); % apply filter to isolate band
dhi = filtfilt(B,A,d); 

% How well does filter work?
plot(t,d);
hold on;
plot(t, dlo, 'r');
plot(t, dhi, 'g');
hold off;
xlabel('Time [s]'); title('Data (blue), Low-freq (red), High-freq (green)');

% Note: issues at edges of data, side-effect of filtering
% center data:

dlo = dlo(N/4: end-N/4-1);
dhi = dhi(N/4: end-N/4-1);
taxis = t(N/4: end-N/4-1);

% Step 2: Compute phase of low freq signal and amplitude envelope of high
% freq signal 

phi = angle(hilbert(dlo)); % phase of low freq signal
amp = abs(hilbert(dhi)); % amplitude env of high freq signal

figure;
subplot(2,1,1); plot(taxis, dlo);
hold on;
plot(taxis, dlo);
hold on;
plot(taxis, phi, 'g');
hold off;
axis tight
xlabel('Time [s]'); title('Low freq and phase');

subplot(2,1,2); plot(taxis, dhi);
hold on;
plot(taxis, amp, 'r');
hold off;
axis tight
xlabel('Time [s]'); title('High freq and amplitude envelope');

% Step 3: Determine if phase and envelope are related. We will do this by
% dividing the phase into binds, find the times where low freq phase lies
% in each bin, and then compute the average amplitude envelope of the high
% freq signal at those times

p_bins = -pi: 0.2: pi;

a_mean = zeros(length(p_bins)-1,1); % vector to hold avg amp env results
p_mean = zeros(length(p_bins)-1,1); % vector to hold center of phase bins


for k = 1:length(p_bins)-1
    pL = p_bins(k); % phase lower limit for this bin
    pR = p_bins(k+1); % phase upper limit for this bin
    indices = find(phi >= pL & phi < pR); % find phase values in this range
    a_mean(k) = mean(amp(indices)); % compute mean amplitude at these phases
    p_mean(k) = mean([pL, pR]); % label the phase bin with the center phase
end

h = max(a_mean) - min(a_mean); % diff bw max and min modulation

% plot the mean envelope vs phase
figure; 
plot(p_mean, a_mean, 'k', 'LineWidth',1);
axis tight
xlabel('Low freq phase'); ylabel('High freq envelope height diff');
title(['Metric h=' num2str(h)]);


%% Conclusion of Demo











































































































