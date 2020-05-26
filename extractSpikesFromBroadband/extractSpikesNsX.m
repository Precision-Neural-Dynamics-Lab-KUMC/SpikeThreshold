% For extracting spikes from .ns5 file saved with Ripple and saving spike waveforms in a .nex file
% Original file: Adam Rouse, 12/4/17
% Modified version: Seng Bum Michael Yoo, 12/07/2017.
% v0.2: Adam Rouse, 12/11/2017
% v0.4: Adam Rouse, 6/30/2019
% v0.4: Adam Rouse, 7/5/2019
% v0.6: Adam Rouse, 7/29/2019

%   Work flow overall: split NSx files via splitNSx.m file, then run by each array.
%   Work flow in current function: make most of things to struct, call functions in dropbox, save back into the HDD.
%   Modification: asof 2017-12-07


function extractSpikesNsX(dataPaths, envInfo, strobeInfo, filtInfo)

% % Path definition (can add, or can direct to that path)
% data_paths.exc_path = fileparts(which( mfilename));
% cd(data_paths.exc_path); cd('../..'); data_paths.ripple_main = pwd;

if filtInfo.use_only_trials && (nargin < 3 || isempty(strobeInfo))
    % Strobe details
    strobeInfo.trial_start_strb    = 'TrialID';  % Special case of TrialID uses all event codes <3000, otherwise works with a single event code for start of trial
    strobeInfo.trial_end_strb      = 6013;  % Right now, the trial_end_code must be immediately before the next trial_start_code, TODO - add more flexibility to time between end and next start code
    strobeInfo.spike_end_strb      = 6013;  % If I put the particular event strobe, then I can get part of the data
    strobeInfo.spike_end_offset    = 0;   % (value: 0 to minus integer) Putting zeros in here will make go through all until spike end code.
end

if nargin < 4 || isempty(filtInfo)
    filtInfo.filt_order = 4;  %4th order filter
    filtInfo.band_limits = [250, 7500]; %bandpass between 250-7500 Hz
    filtInfo.pre_data       = 6;   %Number of data points before trigger time point
    filtInfo.post_data      = 18;  %Number of data points after trigger time point
    filtInfo.peak_window    = 10;  %Number of data points after trigger where waveform peak can occur
    filtInfo.align_spikes   = false;
    filtInfo.throwout_crosstalk = false;
    filtInfo.throwout_large_artifact = false;
end

% if ~isfield(envInfo, 'ch_offset')
%     envInfo.ch_offset = 0;
% end

% path.func_path      = '/home/syoo/Dropbox/HaydenLab_Rouse_Shared/extractSpikesFromBroadband';
% data_paths.input_file_path     = fullfile( data_paths.ripple_main, [envInfo.monkey, '/', envInfo.task, '/' envInfo.date_str] );
% data_paths.median_path   = fullfile( data_paths.file_path, 'SignalQuality' );  %For saving the median values for each channel in a .mat file
% data_paths.save_path     = fullfile( data_paths.file_path, 'Save_data'); 

% cd(data_paths.file_path);


version = '0.4';

fileInfo    = openNSx( [dataPaths.input_file_path , envInfo.ns5_file_name], 'noread');
eventInfo   = openNEV( [dataPaths.input_file_path , envInfo.nev_file_name], 'noread');

%Check to see if electrode channel is actually in data file, remove those that are not, AGR 20200220 
ch_in_data_file = [fileInfo.ElectrodesInfo.ElectrodeID];
for iArr = 1:length(envInfo.channels_to_read_by_array)
    envInfo.channels_to_read_by_array{iArr} = envInfo.channels_to_read_by_array{iArr}(ismember(envInfo.channels_to_read_by_array{iArr}, ch_in_data_file));
    envInfo.file_channels_to_read_by_array{iArr} = find(ismember(ch_in_data_file,envInfo.channels_to_read_by_array{iArr}));
end
envInfo.array_to_fileNum = envInfo.array_to_fileNum(~cellfun(@isempty, envInfo.file_channels_to_read_by_array));
envInfo.channels_to_read_by_array = envInfo.channels_to_read_by_array(~cellfun(@isempty, envInfo.file_channels_to_read_by_array));
envInfo.file_channels_to_read_by_array = envInfo.file_channels_to_read_by_array(~cellfun(@isempty, envInfo.file_channels_to_read_by_array));

if ~isempty(strobeInfo)
%% Event strobe extracting part.
Events = eventInfo.Data.SerialDigitalIO.UnparsedData;
EventTimeSamples    = eventInfo.Data.SerialDigitalIO.TimeStamp;
strobeInfo = find_trial_samples(strobeInfo, Events, EventTimeSamples);
end

numDataPoints = fileInfo.MetaTags.DataPoints(1);

%% Filter construction
%   Sampling frequency (30K)
envInfo.samp_rate = fileInfo.MetaTags.SamplingFreq;
if ~isfield(filtInfo, 'pre_data')
    filtInfo.pre_data       = round(2e-4*envInfo.samp_rate);   %Number of data points before trigger time point
end
if ~isfield(filtInfo, 'post_data')
    filtInfo.post_data      = round(6e-4*envInfo.samp_rate);  %Number of data points after trigger time point
end
if ~isfield(filtInfo, 'peak_window')
    filtInfo.peak_window    = round(3.33e-4*envInfo.samp_rate);   %Number of data points after trigger where waveform peak can occur
end
if ~isfield(filtInfo, 'peak_offset')
    filtInfo.peak_offset = 2; %Number of data points after threshold to put peak of waveform
end
filtInfo.snippet_length = filtInfo.pre_data + filtInfo.post_data + 1;

if length(filtInfo.band_limits) == 1
    [filtInfo.b,filtInfo.a] = butter(filtInfo.filt_order, filtInfo.band_limits(1)/(envInfo.samp_rate/2), 'high');
else
    [filtInfo.b,filtInfo.a] = butter(filtInfo.filt_order, [filtInfo.band_limits(1)/(envInfo.samp_rate/2), filtInfo.band_limits(2)/(envInfo.samp_rate/2)] );
end
if filtInfo.align_spikes && ~isfield(filtInfo, 'interp_factor')
    filtInfo.interp_factor = 4;
end



%% Calculate Median
% cd(data_paths.exc_path);
if ~isfield(filtInfo, 'num_trials_for_median')
filtInfo.num_trials_for_median = 50;
end
calculate_Medians(envInfo, dataPaths, strobeInfo, filtInfo) 


if isempty(strobeInfo)
    trialSize = 1000000;
    trialMargin = 750;
    strobeInfo.trial_start_samp = (0:floor(numDataPoints/trialSize))*trialSize;
    strobeInfo.trial_end_samp = [strobeInfo.trial_start_samp(2:end)+trialMargin, numDataPoints];
    strobeInfo.trial_start_samp(1) = 1;
    strobeInfo.trial_start_samp(2:end) = strobeInfo.trial_start_samp(2:end)-trialMargin;
    strobeInfo.trial_start_samp(end) = max([strobeInfo.trial_start_samp(end),numDataPoints]);   
    strobeInfo.TrialIDs = 1:length(strobeInfo.trial_start_samp);
end

sigQual_name = regexprep(envInfo.ns5_file_name, '.ns\d', '_SigQual.mat');
%% Converting into nex files
%   Create empty nex file in save path and write that at same folder.
% cd( data_paths.save_path );

unique_files = unique(envInfo.array_to_fileNum);



for iFile = unique_files
    curr_ch_count = 0;
nexFileData      = create_blank_nex();
nexFileData.tend = numDataPoints./envInfo.samp_rate;

if ~isempty(eventInfo.Data.SerialDigitalIO.TimeStamp)
    % Events codes are from .nev files.
    nexFileData.markers{1,1}.name = 'Strobed';
    nexFileData.markers{1,1}.varVersion = 100;
    nexFileData.markers{1,1}.values{1,1}.name = 'DIO';
    nexFileData.markers{1,1}.values{1,1}.strings = cellfun(@num2str, num2cell(eventInfo.Data.SerialDigitalIO.UnparsedData), 'Uni', false);
    nexFileData.markers{1,1}.timestamps = eventInfo.Data.SerialDigitalIO.TimeStampSec';
else
    nexFileData = rmfield(nexFileData,{'markers','events'});
end

for iArr = find(envInfo.array_to_fileNum==iFile)  %1:length(envInfo.channels_to_read_by_array)
    spike_snippets  = cell(length(envInfo.channels_to_read_by_array{iArr}),1);
    spike_times     = cell(length(envInfo.channels_to_read_by_array{iArr}),1);
    
    for ch = 1:length(envInfo.channels_to_read_by_array{iArr})
        spike_snippets{ch} = nan(filtInfo.snippet_length,0);
        spike_times{ch} = zeros(0,1,'uint32');
    end
    load([dataPaths.median_path sigQual_name], 'SigQuality')
    curr_index = find(arrayfun(@(x) isequal(x.channels,envInfo.channels_to_read_by_array{iArr}), SigQuality));
    if ~isfield(filtInfo, 'prewhiten_data') ||  ~filtInfo.prewhiten_data
        ChMedians = SigQuality(curr_index).ChMedians;
    else
        ChMedians = SigQuality(curr_index).WhtChMedians;
        WhtMat = SigQuality(curr_index).WhtMat;
    end
    ch_string =  'c:';
    file_channels = sort(envInfo.file_channels_to_read_by_array{iArr},'ascend');
    ch_string = [ch_string, num2str(file_channels(1))];
    colon_flag = true;
    for ch = 2:length(file_channels)
        if (file_channels(ch)-file_channels(ch-1)) > 1
            if ~colon_flag
                ch_string = [ch_string, num2str(file_channels(ch-1))];
            end
            ch_string = [ch_string, ',', num2str(file_channels(ch))];
            colon_flag = true;
        elseif colon_flag
            ch_string = [ch_string, ':'];
            colon_flag = false;
        end
    end
    if ~colon_flag
        ch_string = [ch_string, num2str(file_channels(end))];
    end
    
    for tr = 1:length(strobeInfo.TrialIDs)
        if strobeInfo.trial_end_samp(tr) < numDataPoints
            % Time string
            time_string = ['t:' num2str(strobeInfo.trial_start_samp(tr)) ':' num2str(strobeInfo.trial_end_samp(tr))];
            tempData = openNSx( [dataPaths.input_file_path , envInfo.ns5_file_name], 'read', ch_string, time_string,  'sample');
            
            % Double conversion( % Matlab 2017 specific (?) )
            tempData.Data   = double(tempData.Data);
            tempData        = filter(filtInfo.b, filtInfo.a, tempData.Data');
            if filtInfo.prewhiten_data
                tempData = tempData*WhtMat;
            end
            for ch = 1:length(envInfo.channels_to_read_by_array{iArr})
                if abs(ChMedians(ch)) > 0
                    threshold = -filtInfo.threshold_scale_factor*ChMedians(ch)*0.6745;
                    threshold_crossings = find(tempData(:,ch)<threshold);
                    if exist('trialMargin','var')
                        threshold_crossings = threshold_crossings(threshold_crossings>(trialMargin+1) & threshold_crossings< (size(tempData,1)-trialMargin+1));
                    end
                    if ~isempty(threshold_crossings); threshold_crossings = threshold_crossings([true; (diff(threshold_crossings)>filtInfo.post_data)]); end
                    if ~isempty(threshold_crossings); threshold_crossings = threshold_crossings(threshold_crossings + filtInfo.post_data < size(tempData,1)); end
                    if ~isempty(threshold_crossings); threshold_crossings = threshold_crossings(threshold_crossings - filtInfo.pre_data > 0); end
                    if ~isempty(threshold_crossings)
                        curr_spike_snippets = tempData(threshold_crossings + (-filtInfo.pre_data:filtInfo.post_data), ch);
                        curr_spike_snippets = reshape(curr_spike_snippets,[],filtInfo.snippet_length);
                        
                        if filtInfo.align_spikes
                            peak_times = find_waveform_peaks(curr_spike_snippets', filtInfo.pre_data, filtInfo.peak_window, filtInfo.interp_factor);
                            threshold_crossings = threshold_crossings + floor(peak_times)';
                            threshold_crossings = threshold_crossings(threshold_crossings>(filtInfo.pre_data+filtInfo.peak_offset));
                            threshold_crossings = threshold_crossings((threshold_crossings+(filtInfo.post_data-filtInfo.peak_offset))<size(tempData,1));
                            peak_times_fraction = mod(peak_times,1);
                            curr_spike_snippets = tempData(threshold_crossings + (-(filtInfo.pre_data+filtInfo.peak_offset):(filtInfo.post_data-filtInfo.peak_offset+1)), ch);
                            curr_spike_snippets = reshape(curr_spike_snippets,[],filtInfo.snippet_length+1);
                            for n = 1:size(curr_spike_snippets,1)
                                curr_spike_snippets(n,1:(end-1)) = interp1(1:(filtInfo.snippet_length+1), curr_spike_snippets(n,:), (1:filtInfo.snippet_length)+ peak_times_fraction(n), 'spline');
                            end
                            curr_spike_snippets = curr_spike_snippets(:,1:(end-1));
                        end
                        
                        spike_snippets{ch}  = [spike_snippets{ch}, curr_spike_snippets'];
                        spike_times{ch}     = [spike_times{ch}; strobeInfo.trial_start_samp(tr)+uint32(threshold_crossings)-1];
                    end
                end
            end
        end
        clear tempData
    end
    fprintf('trial: %4.2d   \n', tr );
    
    
    if filtInfo.throwout_crosstalk
        throwout_times = find_crosstalk_times(spike_times, spike_snippets, filtInfo);
        for ch = 1:length(spike_snippets)
            throwout_index = ismember(spike_times{ch}, throwout_times);
            spike_snippets{ch} = spike_snippets{ch}(:, ~throwout_index);
            spike_times{ch} = spike_times{ch}(~throwout_index);
        end
    end
    
    if filtInfo.throwout_large_artifact
        if ~isfield(filtInfo, 'max_number_outlier_spikes')
            filtInfo.max_number_outlier_spikes = 1000;
        end
        if ~isfield(filtInfo, 'outlier_threshold')
            filtInfo.outlier_threshold = 2;
        end
        for ch = 1:length(spike_snippets)
            throwout_index = max(abs(spike_snippets{ch})) > filtInfo.outlier_threshold*prctile( max(abs(spike_snippets{ch})), max([100*(1-filtInfo.max_number_outlier_spikes./size(spike_snippets{ch},2)),90]) );
            spike_snippets{ch} = spike_snippets{ch}(:, ~throwout_index);
            spike_times{ch} = spike_times{ch}(~throwout_index);
        end
    end
    
    
    
    for ch = 1:length(envInfo.channels_to_read_by_array{iArr})
        nexFileData.neurons{curr_ch_count+ch,1}.name = ['sig', num2str(envInfo.channels_to_read_by_array{iArr}(ch), '%03u') , 'U'];
        nexFileData.neurons{curr_ch_count+ch,1}.varVersion = 101;
        nexFileData.neurons{curr_ch_count+ch,1}.wireNumber = curr_ch_count+ch-1;
        nexFileData.neurons{curr_ch_count+ch,1}.unitNumber = 0;
        nexFileData.neurons{curr_ch_count+ch,1}.xPos = 0;
        nexFileData.neurons{curr_ch_count+ch,1}.yPos = 0;
        nexFileData.neurons{curr_ch_count+ch,1}.timestamps = double(spike_times{ch})/envInfo.samp_rate;
        
        nexFileData.waves{curr_ch_count+ch,1}.name = ['sig', num2str(envInfo.channels_to_read_by_array{iArr}(ch), '%03u') , 'U_wf'];
        nexFileData.waves{curr_ch_count+ch,1}.varVersion  = 101;
        nexFileData.waves{curr_ch_count+ch,1}.NPointsWave = size(spike_snippets{ch},1);
        nexFileData.waves{curr_ch_count+ch,1}.WFrequency  = envInfo.samp_rate;
        nexFileData.waves{curr_ch_count+ch,1}.wireNumber  = nexFileData.neurons{curr_ch_count+ch,1}.wireNumber;
        nexFileData.waves{curr_ch_count+ch,1}.unitNumber  = nexFileData.neurons{curr_ch_count+ch,1}.unitNumber;
        nexFileData.waves{curr_ch_count+ch,1}.ADtoMV      = 2.5e-4; %Raw Ripple resolution is 0.25 microvolts or 2.5e-4 mV
        nexFileData.waves{curr_ch_count+ch,1}.MVOffset    = 0;
        nexFileData.waves{curr_ch_count+ch,1}.timestamps = nexFileData.neurons{curr_ch_count+ch,1}.timestamps;
        nexFileData.waves{curr_ch_count+ch,1}.waveforms = nexFileData.waves{curr_ch_count+ch,1}.ADtoMV.*spike_snippets{ch};
    end
    curr_ch_count = length(nexFileData.neurons);
    
end
for k = 1:length(nexFileData.neurons)
    nexFileData.neurons{k}.name((end+1):64) = 0;  %Add null character to name to make 64 character string padded by nulls, important for reading back into plexon offline sorter
    nexFileData.waves{k}.name((end+1):64)   = 0;
end
% writeNex5File(nexFileData, [dataPaths.save_path, envInfo.output_names{iFile} '5']);
writeNexFile(nexFileData, [dataPaths.save_path, envInfo.output_names{iFile}] );
end