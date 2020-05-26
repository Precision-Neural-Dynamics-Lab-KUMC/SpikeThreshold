%Example call script for using extractSpikesFromBroadband package
%Change code with ******* to your path and file information


%*******Root path where files are located
root_file_path = 'C:\Users\Adam\Box\';


%*********Added tools needed: NPMK and HowToReadAndWriteNexAndNex5FilesInMatlab
addpath(genpath([root_file_path 'added_matlab_tools\NPMK-4.5.3.0']))
addpath(genpath([root_file_path 'added_matlab_tools\HowToReadAndWriteNexAndNex5FilesInMatlab']))

%********Add extractSpikesFromBroadBand to path
addpath(genpath('.\extractSpikesFromBroadBand'))

%*****File name to read
envInfo.ns5_file_name	= 'P_20170705_GHIJKLA+A_BB1-64.ns5';

%*****Currently set up assuming .ns5 and .nev have same file name, can also
%be edited here with a name independent of .ns5 name
envInfo.nev_file_name = regexprep( envInfo.ns5_file_name, 'ns\d', 'nev');  


%*********Each element of cell is a group of channels in an array that will have
%their covariance calculated and the common average subtracted away
envInfo.channels_to_read_by_array = {1:16,17:32};

%*******The output .nex file(s) are specified here
%.output_names is the .nex file names to be output
envInfo.output_names = {'P_20170705_GH_out.nex'};
%.array_to_fileNum is which file in .output_names a given array in
%.channels_to_read_by_array should be output 
% All 1's simply outputs all specified channels to the first file name
envInfo.array_to_fileNum = ones(size(envInfo.channels_to_read_by_array));

% envInfo.channels_to_read_by_array = {33:48,49:64};
% envInfo.output_names = {'P_20170705_Iout.nex', 'P_20170705_J_out.nex'};
% envInfo.array_to_fileNum = [1 2 ];



 
root_file_path = 'C:\Users\Adam\Box\';
dataPaths.input_file_path	= [root_file_path 'DataFiles\data_raw\monk_p\20160504_COT_precision\P_CerebusData\'];
dataPaths.input_file_path	= 'C:\DataFiles\data_raw\monk_p\20160504_COT_precision\P_CerebusData\';
dataPaths.median_path		= [root_file_path 'DataFiles\data_processed\monk_p\20160504_COT_precision\SignalQuality\'];
dataPaths.save_path			= [root_file_path 'DataFiles\data_processed\monk_p\20160504_COT_precision\BB_to_Spikes\'];



filtInfo.overwrite_median = true;
filtInfo.num_trials_for_median = 50;
filtInfo.filt_order = 4;  % 4th order filter
filtInfo.band_limits = [250, 5000]; % bandpass between 250-7500 Hz
filtInfo.pre_data       = 12;    % Number of data points before trigger time point
filtInfo.post_data      = 38;   % Number of data points after trigger time point
filtInfo.align_spikes   = false; %Align spikes based on peak to get better alignment of waveforms
filtInfo.interp_factor  = 4;    %Align spikes based on peak to get better alignment of waveforms
filtInfo.peak_window    = 10;   %Number of data points after trigger where waveform peak can occur
filtInfo.peak_offset    = 2;    %Number of data points after threshold to put peak of waveform
filtInfo.prewhiten_data = true;
filtInfo.throwout_crosstalk = true;
filtInfo.num_ct_channels = 5;
filtInfo.std_per_chan   = 5;
filtInfo.throwout_large_artifact = true;
filtInfo.max_number_outlier_spikes = 1000;
filtInfo.outlier_threshold = 2;
filtInfo.threshold_scale_factor = 10;  %Threshold scale factor

filtInfo.use_only_trials = false;
filtInfo.median_window = 1000;  %1000 ms window for median (only used when not using actual trials)

% strobeInfo.trial_start_strb    = 'TrialID';  % Special case of TrialID uses all event codes <3000, otherwise works with a single event code for start of trial
% strobeInfo.trial_end_strb      = 6013;  % Right now, the trial_end_code must be immediately before the next trial_start_code, TODO - add more flexibility to time between end and next start code
% strobeInfo.spike_end_strb      = 6013;  % If I put the particular event strobe, then I can get part of the data
% strobeInfo.spike_end_offset    = 0;   % (value: 0 to minus integer) -integer is number of samples before spike_end_strb, putting 0 here will go all the until spike_end_strb.

% strobeInfo.trial_start_strb    = 'TrialID';  % Special case of TrialID uses all event codes <3000, otherwise works with a single event code for start of trial
% strobeInfo.trial_end_strb      = 6253;  % Right now, the trial_end_code must be immediately before the next trial_start_code, TODO - add more flexibility to time between end and next start code
% strobeInfo.spike_end_strb      = 6120;  % If I put the particular event strobe, then I can get part of the data
% strobeInfo.spike_end_offset    = -10;   % (value: 0 to minus integer) -integer is number of samples before spike_end_strb, putting 0 here will go all the until spike_end_strb.

strobeInfo = [];

extractSpikesNsX(dataPaths, envInfo, strobeInfo, filtInfo)
