function filtInfo = defaultFiltInfo()
%Create default filter info for extractSpikesREC
filtInfo.overwrite_median = true;
filtInfo.num_trials_for_median = 50;
filtInfo.filt_order = 4;  % 4th order filter
filtInfo.band_limits = [250, 5000]; % bandpass between 250-7500 Hz
filtInfo.time_pre       = 175;    % Amount of time before trigger for snippet (microseconds)
filtInfo.time_post      = 625;    % Amount of time after trigger for snippet (microseconds)
filtInfo.time_peak_excl = 625;   %Minimum time from previous threshold crossing that the next spike can occur
filtInfo.time_req_baseline = 175;  %Minimum time signal must be below threshold crossing before next spike can occur
filtInfo.align_spikes   = true; %Align spikes based on peak to get better alignment of waveforms
filtInfo.interp_factor  = 4;    %Interpolation factor for alignment of spikes 
filtInfo.peak_time      = 150;   %Time after trigger where waveform peak can occur  (microseconds)
filtInfo.peak_offset    = 2;    %Number of data points after threshold to put peak of waveform
filtInfo.prewhiten_data = true;
filtInfo.throwout_crosstalk = true;
filtInfo.num_ct_channels = 5;
filtInfo.std_per_chan   = 5;
filtInfo.throwout_large_artifact = true;
filtInfo.max_number_outlier_spikes = 1000;
filtInfo.outlier_threshold = 2;
filtInfo.threshold_scale_factor = 7.5;  %Threshold scale factor

filtInfo.use_only_trials = false;
filtInfo.median_window = 1000;  %1000 ms window for median (only used when not using actual trials)
end

