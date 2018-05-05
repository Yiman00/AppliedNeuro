% Create some data
dt = 0.001; % time step [s]
t = (dt:dt:3); % duration [s]
T = max(t); % max time
N = length(t); % total number of indices
x = sin(2.0*pi*t*10) + randn(1,N); % signal + noise data

% Fourier transform, decomposes signal into sinusoids 
% Nyquist frequency (sampling frequency/2), in steps of 1/T
% - Nyquist frequency: time step: + Nyquist frequency
fj = (-1/dt/2: 1/T :1/dt/2-1/T); 

% For each freq, compute the Fourier transform and save as X
X = zeros(1, length(fj)); % pre-allocate
for j = 1:length(fj)
    X(j) = sum( x.* exp(-2*pi*1i*fj(j)*t) ); 
end

% Compute the power spectrum
pow = 2*dt^2/T * X.*conj(X); 
subplot(3,1,1); 
plot(fj, 10*log10(real(pow))); axis tight; xlim([0 20]); ylim([-50 0]); 
xlabel('Freq Hz'); ylabel('Power'); % [V^2], spectral density (V^2/ Hz)

% Plot with built-in FFT function
X = fft(x); % y-axis

% make x-axis for built-in FFT plot (freq)
df = 1/T;
fNQ = 1/dt/2; % Nyquist freq
faxis = (-fNQ : df : fNQ-df); 

% Re-compute the power spectrum
pow = 2*dt^2/T * X.*conj(X); 

% Plot on decibel scale
subplot(3,1,2);
plot(faxis, 10*log10(fftshift(pow))); axis tight; xlim([0 20]); ylim([-50 0]); 
xlabel('Freq Hz'); ylabel('Power'); 

% Plot power spectrum
subplot(3,1,3);
periodogram(x, [], N, 1/dt); xlim([ 0 20]); ylim([-50 0]); 

% V(X), V(Y) => V(X,Y) co-variance
% Compute co-variance 

% Make random noisy signals
dt = 0.001; % sampling interval
t = dt:dt:1; % [s]
x = 0.2*randn(1, length(t));
y = 0.2*randn(1, length(t)); 

% plot two signals
figure; 
plot(t, x); 
hold on;
plot(t, y, 'g');
hold off;
xlabel('Time [s]'); 

















    
    
    
    
    

