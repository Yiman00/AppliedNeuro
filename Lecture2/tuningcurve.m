% Tuning Curve

ntrials = numel(tonespike); % number of trials
Feature = NaN(ntrials,1); % pre-allocation 
for ii = 1:ntrials 
    values = tonespike(ii).stimvalues; 
    Feature(ii) = values(1); % (Hz)
end

sel = Feature == 250; % select all trials with tone frequency 
spiketime = [tonespike(sel).spiketime]; 
nrep = sum(sel); % number of stimulus repeats

onset = 300; % tone onset (ms)
offset = 450; % tone offset (ms);
sel = spiketime>=onset & spiketime<=offset; % spikes during stimulus 
Nspikes = sum(sel); % number of spikes

duration = offset - onset; % ms
Fs = 1000; % sampling rate
firing_rate = Nspikes/duration*Fs/nrep; % in spikes/s or Hz


