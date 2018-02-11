% Tuning Curve for Freq

uFeature = unique(Feature); % the unique frequencies presented during the recordings
nFeature = numel(uFeature); % the number of freqs used
firingRate = NaN(nFeature,1); % pre-allocation

Fs = 1000; % sampling rate (Hz) 
onset = 300; 
offset = 450;
duration = offset-onset;

for ii = 1:nFeature
    sel = Feature == uFeature(ii); % select all trials with tone freq of
    % uFeature(ii)
    spiketime = [tonespike(sel).spiketime];
    nrep = sum(sel); % number of stimulus repeats
    sel = spiketime>=onset & spiketime<=offset; % stimulus in duration
    Nspikes = sum(sel); % number of spikes
    firingRate(ii) = Nspikes/duration*Fs/nrep; % spikes/s
end

% Plot Data
h = semilogx(uFeature,firingRate); % Plot
set(h, 'Color','k','LineWidth',2,'LineStyle','-'); 
set(gca, 'XTick', uFeature(1:2:end), 'XTickLabel', round(uFeature(1:2:end))); box off; 
xlabel('Frequency (Hz)');
ylabel('Firing Rate (Hz)');

axis square;
xlim([0.9*min(uFeature) 1.1*max(uFeature)]);
ylim([0 60]);










    