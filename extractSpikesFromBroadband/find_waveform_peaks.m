function peak_times = find_waveform_peaks(spike_waveforms, pre_data, peak_window, interp_factor)


%Preallocate
%     interp_spike = zeros((peak_window+1)*interp_factor, size(spike_waveforms,2));
%     %Interpolate around the alignment_index to better estimate the exact time
%     %point where the peak or trough occurred,  interpolation is done
%     %interp_offset samples before and after alignment_index, with
%     %interp_factor more samples
%     for n = 1:size(spike_waveforms,2)
%        interp_spike(:,n) = interp(spike_waveforms((pre_data+1):(pre_data+peak_window+1),n), interp_factor);
%     end
    interp_spike = interp1(1:(peak_window+1), spike_waveforms((pre_data+1):(pre_data+peak_window+1),:), 1:(1/interp_factor):(peak_window+1), 'spline');
    
    [~, peak_times] = min(interp_spike,[],1);
 peak_times = peak_times./interp_factor;

