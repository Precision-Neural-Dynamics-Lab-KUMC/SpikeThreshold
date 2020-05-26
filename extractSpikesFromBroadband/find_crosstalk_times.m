function throwout_times = find_crosstalk_times(spike_times, spike_snippets, filtInfo)
%% fine_noise_times - Identify noise traces based on number spikes happening at same time across channels and the magnitude of the traces
% Adam Rouse, 12/5/17
% num_ct_channels - Number of  channels with signals at same time required to throwout
% std_per_chan         - multiple of std of all waveforms for setting threshold
if ~isfield(filtInfo, 'num_ct_channels')
    filtInfo.num_ct_channels = 5;
end
if ~isfield(filtInfo, 'std_per_chan')
    filtInfo.std_per_chan = 5;
end


all_spike_times = cat(1,spike_times{:});
spike_amplitudes = cellfun(@(x) max(abs(x)), spike_snippets, 'Uni', false);
all_spike_amplitudes = cat(2,spike_amplitudes{:})';

[all_spike_times, sort_order] = sort(all_spike_times);
all_spike_amplitudes = all_spike_amplitudes( sort_order );

unique_times = unique(all_spike_times);

spike_counts = histcounts(all_spike_times,[unique_times;Inf]);

ci_times_i = find(spike_counts>=filtInfo.num_ct_channels);

sum_spike_amplitudes = zeros(length(ci_times_i),1);

for k = 1:length(ci_times_i)
    sum_spike_amplitudes(k) = sum(all_spike_amplitudes(all_spike_times==unique_times(ci_times_i(k))));
end

threshold = filtInfo.std_per_chan*std(all_spike_amplitudes)*filtInfo.num_ct_channels;

throwout_times = unique_times(ci_times_i(sum_spike_amplitudes>threshold));


% for ch = 1:length(spike_snippets)
%     good_spike_snippets{ch} = spike_snippets{ch}(:, ~ismember(spike_times{ch}, throwout_times));
%     bad_spike_snippets{ch} = spike_snippets{ch}(:, ismember(spike_times{ch}, throwout_times));
% end



% figure
% plot(bad_spike_snippets{12}(:,1:10:end), 'Color', [0.5 0.5 0.5])
% hold on
% plot(good_spike_snippets{12}(:,1:10:end), 'Color', [0 0 0])




