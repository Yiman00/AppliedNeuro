% E M Izhikevich: Network Model of Mammalian Cortex

close all
clc
clear 

% INPUT PARAMETERS   =====================================================
   Ne = 800; Ni = 200;         % number of excitatory & inhibitory neurons
   Nt = 1000;                  % number of time steps
   numE = 10; numI = Ne + 10;  % indices for one excitatory and one inhibitory neuron
   
% Model parameters   =====================================================
   re = rand(Ne,1); ri = rand(Ni,1);     % random numbers
   a = [0.02 * ones(Ne,1); 0.02 + 0.08*ri];
   b = [0.20   *ones(Ne,1); 0.25 - 0.05*ri];
   c = [-65+15*re.^2; -65*ones(Ni,1)];
   d = [8-6*re.^2   ; 2*ones(Ni,1)];
   S = [0.53*rand(Ne+Ni,Ne), -rand(Ne+Ni,Ni)];   % coupling strengths

   v = -65*ones(Ne+Ni,1); % Initial values of v
   u =  b.*v;             % Initial values of u
   firings = [];          % spike timings

   vE = zeros(Nt,1); vI = zeros(Nt,1);   % membrane potential of 2 neurons

   
% Time Evolution of Systems  ==============================================
for t = 1:Nt 
   I = [5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
   fired = find(v>=30);               % indices of spikes
   firings = [firings; t+0*fired,fired];  % time steps  / fired neurons
   v(fired) = c(fired);               % membrane potential
   u(fired) = u(fired)+d(fired);      % recovery potential
   I = I+sum(S(:,fired),2);           % thalamic + synaptic input
   
   v = v+0.5*(0.04*v.^2+5*v+140-u+I); % HALF-STEP: step 0.5 ms
   v = v+0.5*(0.04*v.^2+5*v+140-u+I); %  for numerical stability
   u = u+a.*(b.*v-u); 
   
  vE(t) = v(numE);                    % excitatory neuron 
  vI(t) = v(numI);                    % inhibitory neuron
end

   vE(vE > 30) = 30;
   vI(vI > 30) = 30;


tS = 1:Nt; nF = zeros(Nt,1);        % fired neurons at each time step
for t = 1 : Nt
   nF(t) = sum((firings(:,1) == t));
end

% GRAPHICS ===============================================================

figure(1)    % raster plot of firing neurons / % firing neurons
   set(gcf,'units','normalized','Position',[0.1 0.4 0.32,0.45]);
   set(gca,'fontsize',12);
   col = [1 0 0];
   hP = subplot(2,1,1);   % raster plot
   plot(firings(:,1),firings(:,2),'.','color',col);
   set(hP, 'Position',[0.1300 0.52 0.7750 0.4]);
   ylabel('neuron #')
   
   subplot(2,1,2);       % percentage of firong neurons each time step
   plot(tS,nF*100/(Ne+Ni),'color',col);
   xlabel('time steps  [ms]');
   ylabel('% firing spikes')
   grid on; box on;
   
   
figure(2)   % membrane potential: excitatory neuron & inhibitory neuron   
  set(gcf,'units','normalized','Position',[0.3 0.1 0.32,0.30]);
  set(gca,'fontsize',12);
  plot(1:Nt,vE,'b');
  hold on
  plot(1:Nt,vI,'r');
  xlabel('time steps  [ms]');
  ylabel('Membrane Voltage  [mV]')
  grid on; box on;
  legend('excitatory','inhibitory'); 