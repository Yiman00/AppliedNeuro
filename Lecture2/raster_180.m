% Raster plot of the 180th trial

t = tonespike(180).spiketime; % Spike timings in the 180th trial
nspikes = numel(t); % number of spikes 
for ii = 1:nspikes % for every spike
    line([t(ii) t(ii)], [179.5 180.5], 'Color', 'k'); % draw a black
    % vertical line at time t(x) and at trial 180 (y) 
end

xlabel('Time (ms)'); % Time in milliseconds 
ylabel('Trial number'); 