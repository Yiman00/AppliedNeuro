% spike002.m

% 161012
% Ian Cooper
% School of Physics, University of Sydney

% E M Izhikevich: Simple Models of Spiking Neurons

%close all
clc
clear all

% INPUT PARAMETERS   =====================================================

flagS = 1;    % select spike train pattern
N = 5000;     % number of grid points     
tMax = 100;   % max time [ms]

% model parameters a b c d I
coeffM = zeros(20,5);
coeffM(1,:)  = [0.02      0.4     -65      2       10  ] ;
%coeffM(1,:)  = [0.02      0.20     -65      6       20  ] ;  % tonic spiking
coeffM(2,:)  = [0.02      0.25    -65      6       1] ; % phasic spiking
coeffM(3,:)  = [0.02      0.20     -50      2       15  ] ;  % tonic bursting
coeffM(4,:)  = [0.02      0.25    -55     0.05     0.6 ] ;  % phasic bursting
coeffM(5,:)  = [0.02      0.2     -55     4        10  ] ;  % mixed mode
coeffM(6,:)  = [ 0.01      0.2     -65     8        30 ] ;  % spike frequency adaptation 
coeffM(7,:)  = [0.02      -0.1    -55     6        0   ] ;  % Class 1
coeffM(8,:)  = [0.2       0.26    -65     0        0   ] ;  % Class 2
coeffM(9,:)  = [0.02      0.2     -65     6        7   ] ;  % spike latency
coeffM(10,:) = [0.05      0.26    -60     0        0   ] ;  % subthreshold oscillations
coeffM(11,:) = [0.1       0.26    -60     -1       0   ] ;  % resonator
coeffM(12,:) = [ 0.02      -0.1    -55     6        0  ] ;  % integrator
coeffM(13,:) = [0.03      0.25    -60     4        0   ] ;  % rebound spike
coeffM(14,:) = [0.03      0.25    -52     0        0   ] ;  % rebound burst
coeffM(15,:) = [0.03      0.25    -60     4        0   ] ;  % threshold variability
coeffM(16,:) = [1         1.5     -60     0      -65   ] ;  % bistability
coeffM(17,:) = [ 1       0.2     -60     -21      0    ] ;  % DAP
coeffM(18,:) = [ 0.02      1       -55     4        0  ] ;  % accomodation
coeffM(19,:) = [ -0.02      -1      -60     8       80 ] ;  % inhibition-induced spiking
coeffM(20,:) = [ -0.026     -1      -45     0        80] ;  % inhibition-induced bursting
      
   
   a = coeffM(flagS,1);
   b = coeffM(flagS,2);
   c = coeffM(flagS,3);
   d = coeffM(flagS,4);
   I = coeffM(flagS,5);
    
% Plot Title    
   tm = [' 1: Tonic Spiking                '; ...
         ' 2: Phasic Spiking               '; ...
         ' 3: Tonic Bursting               '; ...
         ' 4: Phasic Bursting              '; ...
         ' 5: Mixed Mode                   '; ...
         ' 6: Spike Frequency Adaptation   '; ...
         ' 7: Class 1                      '; ...
         ' 8: Class 2                      '; ...
         ' 9: Spike Latency                '; ...
         '10: Subthreshold Oscillations    '; ...
         '11: Resonator                    '; ...
         '12: Integrator                   '; ...
         '13: Rebound Spike                '; ...
         '14: Rebound Burst                '; ...
         '15: Threshold Variability        '; ...
         '16: Bistability                  '; ...
         '17: DAP                          '; ...
         '18: Accomodation                 '; ...
         '19: Inhibition-induced Spiking   '; ...
         '20: Inhibition-induced Bursting  '; ];
   
   
% SOLVE EQUATIONS   ======================================================   
       
cc(1) = 0.04; cc(2) = 5; cc(3) = 140;     % D.E. coefficients

t = linspace(0,tMax,N);                   % time steps [ms]
v = zeros(N,1); vM = zeros(N,1);          
u = zeros(N,1); uM = zeros(N,1);          % recovery variable
Iext = zeros(N,1);

Iext(N/10:end) = I;
dt = t(2)-t(1);

v(1) = -65;
vM(1) = -65;
u(1) = b * v(1);
uM(1) = u(1);


for n = 1 : N-1
  v(n+1) = v(n) + dt*( cc(1)*v(n)^2 + cc(2)*v(n) + cc(3) + Iext(n) - u(n) );
  
  u(n+1) = u(n) + dt*a*( b*v(n) - u(n) );
  
  if v(n+1) > 30;
     vM(n+1) = 30;
     v(n+1) = c;
     u(n+1) = u(n+1) + d;
     
  else
     vM(n+1) = v(n+1);
     uM(n+1) = u(n+1);
  end
  
   
  
end
figure(2)
set(gcf,'units','normalized','Position',[0.3 0.32 0.32,0.38]);
set(gca,'fontsize',12);
hSubplot = subplot(2,1,1);
set(hSubplot,'position',[0.13 0.46 0.78 0.45]);
xP = t; yP = vM;
plot(xP,yP,'b','lineWidth',2);     
ylabel('v_M  [mV]');
% title(tm(flagS,:));
grid on
set(gca,'fontsize',12);


hSubplot = subplot(2,1,2);
set(hSubplot,'position',[0.13 0.15 0.78 0.20]);
xP = t; yP = Iext;
plot(xP,yP,'r','lineWidth',2);
xlabel('time t  [ms]');
ylabel('c_5 I_{ext}  [mV]');
axis([0 tMax 0 1.2*max(Iext)]);
grid on
set(gca,'fontsize',12);

    