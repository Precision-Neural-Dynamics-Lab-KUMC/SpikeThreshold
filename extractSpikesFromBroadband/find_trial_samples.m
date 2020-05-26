function strobeInfo = find_trial_samples(strobeInfo, Events, EventTimeSamples)
%%find_trial_samples - find indexes for when trials start and end to extract spikes
%Adam Rouse, 12/11/2017
EventItemCounts     = length(Events);
if strcmpi(strobeInfo.trial_start_strb, 'TrialID')
    strobeInfo.trial_start_idx = find(Events<3000);
else
    strobeInfo.trial_start_idx = find(Events==strobeInfo.trial_start_strb);
end
strobeInfo.trial_end_idx    = [strobeInfo.trial_start_idx(2:end)-1; EventItemCounts];
strobeInfo.trial_start_idx  = strobeInfo.trial_start_idx(Events(strobeInfo.trial_end_idx)==strobeInfo.trial_end_strb);
strobeInfo.trial_end_idx    = strobeInfo.trial_end_idx(Events(strobeInfo.trial_end_idx)==strobeInfo.trial_end_strb);
if strcmpi(strobeInfo.trial_start_strb, 'TrialID')
    strobeInfo.TrialIDs = Events(strobeInfo.trial_start_idx);
end
strobeInfo.trial_start_samp = EventTimeSamples(strobeInfo.trial_start_idx);

% Various matching methods is possible.
%   This version is asof 2017-12-07
if strobeInfo.trial_end_strb == strobeInfo.spike_end_strb
    strobeInfo.spike_end_idx = strobeInfo.trial_end_idx;
else
    for tr = 1:length(strobeInfo.trial_start_idx)
        curr_spike_end_idx = strobeInfo.trial_start_idx(tr) - 1 + find(Events(strobeInfo.trial_start_idx(tr):strobeInfo.trial_end_idx(tr))==strobeInfo.spike_end_strb, 1, 'last');
        if ~isempty(curr_spike_end_idx)
            strobeInfo.spike_end_idx(tr) = curr_spike_end_idx;
        else
            strobeInfo.spike_end_idx(tr) = nan;
        end
    end
    strobeInfo.trial_start_idx = strobeInfo.trial_start_idx(~isnan(strobeInfo.spike_end_idx));
     strobeInfo.trial_end_idx = strobeInfo.trial_end_idx(~isnan(strobeInfo.spike_end_idx));
     strobeInfo.TrialIDs = strobeInfo.TrialIDs(~isnan(strobeInfo.spike_end_idx));
     strobeInfo.spike_end_idx = strobeInfo.spike_end_idx(~isnan(strobeInfo.spike_end_idx));
     
end
if length(strobeInfo.spike_end_idx) ~= length(strobeInfo.trial_end_idx)
    error('not same number of spike_end_code events as trial events');    
end
strobeInfo.trial_end_samp = EventTimeSamples(strobeInfo.spike_end_idx) + strobeInfo.spike_end_offset;