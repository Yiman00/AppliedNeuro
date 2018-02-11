% Raster Plot for all trials

ntrials = numel(tonespike); % number of trials 
for jj = 1:ntrials
    t = tonespike(jj).spiketime; % spike timings in the jjth trial
    nspikes = numel(t); % number of spikes
    for ii = 1:nspikes % for every spike
        line([t(ii) t(ii)], [jj-0.5 jj+0.5], 'Color','k');
        % draw a black vertical line of length 1 at time t(x) and at trial 
        % jj (y)
    end
end

xlabel('Time (ms)');
ylabel('Trial number'); 

