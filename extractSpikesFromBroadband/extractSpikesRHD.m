% For extracting spikes from .ns5 file saved with Ripple and saving spike waveforms in a .nex file
% Original file: Adam Rouse, 12/4/17
% Modified version: Seng Bum Michael Yoo, 12/07/2017.
% v0.2: Adam Rouse, 12/11/2017
% v0.4: Adam Rouse, 6/30/2019
% v0.4: Adam Rouse, 7/5/2019
% v0.6: Adam Rouse, 7/29/2019
% v0.7: Adam Rouse, 3/20/2021
% v0.8: Adam Rouse, 7/01/2021
% v0.9: Adam Rouse, 4/15/2022
% v1.0: Adam Rouse, 2/1/2023
% v1.1: Willy Lee, 6/1/2023 adapt for Intan systen
% v1.2: Adam Rouse, 8/13/2024, changed file opening to more explicitly use port and channel number rather than relying on dir
% v1.3: Xavier Scherschligt implemented Dr. Rouse's EventDropTest_final.m (L:58-273),  9/10/24
% v2.0: Adam Rouse, 5/9/2025 Merging back to a single extractSpikes call


function extractSpikesRHD(dataPaths, envInfo, dataBlocks, filtInfo)

if isempty(which('read_Intan_RHD2000_file.m'))
    addpath('./IntanFileReading/')
end

data_strut = readRHDExtractedDataFile([dataPaths.input_file_path 'amp-A-000.dat']);



%Check to see if electrode channel is actually in data file, remove those that are not, AGR 20200220
%add port and signal band determination to check channel numbers in Intan sys, need to hard code the ports being used, WL 20230531
file_names =  dir([dataPaths.input_file_path, '*.dat']);
ch_in_data_file =[];
for k = 1:length(file_names)
    tmp = str2num(file_names(k).name((regexpi(  file_names(k).name, '\d*')):(regexpi(  file_names(k).name, '.dat')-1)));
    num_port = file_names(k).name((regexpi(  file_names(k).name, '\d*')-2));

    if ~isempty(tmp)
        sig_band = file_names(k).name(1:3);
        if length(file_names) < 384
            %when recording 64 ch ports
            if all(sig_band == 'amp')
                if num_port == 'A'
                    ch_in_data_file(end+1) = tmp+1;
                elseif num_port == 'B'
                    ch_in_data_file(end+1) = tmp+33;
                elseif num_port == 'C'
                    ch_in_data_file(end+1) = tmp+65;
                elseif num_port == 'D'
                    ch_in_data_file(end+1) = tmp+97;
                end
            end
        else
            %when recording 128 ch ports
            if all(sig_band == 'amp')
                if num_port == 'A'
                    ch_in_data_file(end+1) = tmp+1;
                elseif num_port == 'B'
                    ch_in_data_file(end+1) = tmp+129;
                elseif num_port == 'C'
                    ch_in_data_file(end+1) = tmp+257;
                end
            end
        end
    end
end
ch_in_data_file = sort(ch_in_data_file);
% ch_in_data_file = [fileInfo.ElectrodesInfo.ElectrodeID];
for iArr = 1:length(envInfo.channels_to_read_by_array)
    envInfo.channels_to_read_by_array{iArr} = envInfo.channels_to_read_by_array{iArr}(ismember(envInfo.channels_to_read_by_array{iArr}, ch_in_data_file));
    envInfo.file_channels_to_read_by_array{iArr} = find(ismember(ch_in_data_file,envInfo.channels_to_read_by_array{iArr}));
end
envInfo.array_to_fileNum = envInfo.array_to_fileNum(~cellfun(@isempty, envInfo.file_channels_to_read_by_array));
envInfo.channels_to_read_by_array = envInfo.channels_to_read_by_array(~cellfun(@isempty, envInfo.file_channels_to_read_by_array));
envInfo.file_channels_to_read_by_array = envInfo.file_channels_to_read_by_array(~cellfun(@isempty, envInfo.file_channels_to_read_by_array));



numDataPoints = length(data_strut.fields.data);


%% Filter construction
disp("In extractSpikesRHD - Constructing Filter")
%   Sampling frequency (30K)
if ~isfield(envInfo, 'samp_rate')
    envInfo.samp_rate = data_strut.clockrate;
elseif envInfo.samp_rate ~= data_strut.clockrate
    error('Neural data and events have different sampling rate!')
end
if ~isfield(filtInfo, 'pre_data')
    if isfield(filtInfo, 'time_pre')
        filtInfo.pre_data       = ceil(filtInfo.time_pre*1e-6*envInfo.samp_rate);   %Number of data points before trigger time point
    else
        filtInfo.pre_data       = ceil(1.75e-4*envInfo.samp_rate);   %Number of data points before trigger time point
    end
end
if ~isfield(filtInfo, 'post_data')
    if isfield(filtInfo, 'time_post')
        filtInfo.post_data      = ceil(filtInfo.time_post*1e-6*envInfo.samp_rate);   %Number of data points before trigger time point
    else
        filtInfo.post_data      = ceil(6.25e-4*envInfo.samp_rate);  %Number of data points after trigger time point
    end
end
if ~isfield(filtInfo, 'peak_excl_data')
    if isfield(filtInfo, 'time_peak_excl')
        filtInfo.peak_excl_data     = ceil(filtInfo.time_peak_excl*1e-6*envInfo.samp_rate);   %Number of data points from previous threshold crossing before next spike threshold crossing can occur
    else
        filtInfo.peak_excl_data     = filtInfo.post_data ;  %Number of data points from previous threshold crossing before next spike threshold crossing can occur
    end
end
if ~isfield(filtInfo, 'req_base_data')
    if isfield(filtInfo, 'time_req_baseline')
        filtInfo.req_base_data     = ceil(filtInfo.time_req_baseline*1e-6*envInfo.samp_rate);   %Number of data points of baseline (not past threshold) needed before next spike threshold crossing can occur
    else
        filtInfo.req_base_data     = filtInfo.pre_data ;   %Number of data points of baseline (not past threshold) needed before next spike threshold crossing can occur
    end
end

if ~isfield(filtInfo, 'peak_window')
    if isfield(filtInfo, 'peak_time')
        filtInfo.peak_window    = ceil(filtInfo.peak_time*1e-6*envInfo.samp_rate);   %Number of data points after trigger where waveform peak can occur
    else
        filtInfo.peak_window    = ceil(3e-4*envInfo.samp_rate);   %Number of data points after trigger where waveform peak can occur
    end
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
disp("In extractSpikesRHD - Calculating Median")
% cd(data_paths.exc_path);
if ~isfield(filtInfo, 'num_trials_for_median')
    filtInfo.num_trials_for_median = 50;
end
calculate_MediansRHD(envInfo, dataPaths, [], filtInfo)



if isempty(dataBlocks)
    trialSize = 10000000;
    trialMargin = 750;
    dataBlocks.block_start_samp = (0:floor(numDataPoints/trialSize))*trialSize;
    dataBlocks.block_end_samp = [dataBlocks.block_start_samp(2:end)+trialMargin, numDataPoints];
    dataBlocks.block_start_samp(1) = 1;
    dataBlocks.block_start_samp(2:end) = dataBlocks.block_start_samp(2:end)-trialMargin;
    dataBlocks.BlockIDs = 1:length(dataBlocks.block_start_samp);
end


sigQual_name = [envInfo.rec_file_name, '_SigQual.mat'];
%% Converting into nex files
disp("In extractSpikesRHD - Converting into nex files")
%   Create empty nex file in save path and write that at same folder.
% cd( data_paths.save_path );

unique_files = unique(envInfo.array_to_fileNum);


for iFile = unique_files
    curr_ch_count = 0;
    nexFileData      = create_blank_nex();
    nexFileData.tend = numDataPoints./envInfo.samp_rate;

    if isfield(envInfo, 'EventTimeSamples') &&  ~isempty(envInfo.Events) && ~isempty(envInfo.EventTimeSamples)
       
        nexFileData.markers{1,1}.name = 'Strobed';
        nexFileData.markers{1,1}.varVersion = 100;
        nexFileData.markers{1,1}.values{1,1}.name = 'DIO';
        nexFileData.markers{1,1}.values{1,1}.strings = cellfun(@num2str, num2cell(envInfo.Events), 'Uni', false);
        nexFileData.markers{1,1}.timestamps = double(envInfo.EventTimeSamples)./double(envInfo.samp_rate);
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



        %Loops through time by dataBlocks using start and stop from
        %dataBlocks.block_start_samp & dataBlocks.block_end_samp
        %Not actual behavioral trials, just consecutive chunks of data
        for tr = 1:length(dataBlocks.BlockIDs)
            if dataBlocks.block_end_samp(tr) <= numDataPoints
                for ch = 1:length(SigQuality(iArr).channels)
                    if SigQuality(curr_index).port(ch) == 'A'
                        ch_offset = 0;
                    elseif SigQuality(curr_index).port(ch) == 'B'
                        ch_offset = 128;
                    elseif SigQuality(curr_index).port(ch) == 'C'
                        ch_offset = 256;
                    end
                    curr_file_index = find(arrayfun(@(x) strcmpi(x.name, ['amp-' SigQuality(curr_index).port(ch) '-' num2str(SigQuality(curr_index).channels(ch)-ch_offset-1,'%03i') '.dat']), file_names));

                    fileData = readRHDExtractedDataFile([file_names(curr_file_index).folder, '\', file_names(curr_file_index).name],dataBlocks.block_start_samp(tr),dataBlocks.block_end_samp(tr));

                    if ch == 1
                        tempData = zeros(size(fileData.fields.data,1),32,'double');
                    end
                    tempData(:,ch) = double(fileData.fields.data);
                    
                    clear fileData
                end
                tempData        = filter(filtInfo.b, filtInfo.a, tempData);
                if filtInfo.prewhiten_data
                    tempData = tempData*WhtMat;
                end
                for ch = 1:length(envInfo.channels_to_read_by_array{iArr})
                    if abs(ChMedians(ch)) > 0
                        threshold = -filtInfo.threshold_scale_factor*ChMedians(ch)*0.6745;
                        threshold_crossings = find(tempData(:,ch)<threshold);
                        if exist('trialMargin','var')
                            if tr == 1  %if first chunk of data, no trialMargin at front
                                threshold_crossings = threshold_crossings(threshold_crossings< (size(tempData,1)-trialMargin+1));
                            elseif tr == length(dataBlocks.BlockIDs) %if last chunk of data, no trialMargin at end
                                threshold_crossings = threshold_crossings(threshold_crossings>(trialMargin+1));
                            else %All other chunks, remove margin at front and back
                                threshold_crossings = threshold_crossings(threshold_crossings>(trialMargin+1) & threshold_crossings< (size(tempData,1)-trialMargin+1));
                            end
                        end
                        %Remove any threshold crossings within the window exclusion
                        if ~isempty(threshold_crossings); threshold_crossings = threshold_crossings([true; (diff(threshold_crossings)>filtInfo.req_base_data)]); end
                        %Remove any threshold crossings within the window exclusion
                        if length(threshold_crossings)>2
                            short_crossings = [false; diff(threshold_crossings)<filtInfo.peak_excl_data & [true; diff(threshold_crossings(1:(end-1)))>=filtInfo.peak_excl_data]];
                            while any(short_crossings)
                                threshold_crossings = threshold_crossings(~short_crossings);
                                short_crossings = [false; diff(threshold_crossings)<filtInfo.peak_excl_data & [true; diff(threshold_crossings(1:(end-1)))>=filtInfo.peak_excl_data]];
                            end
                        end
                        %Remove any snippets that don't have data at the end of the file
                        if ~isempty(threshold_crossings); threshold_crossings = threshold_crossings(threshold_crossings + filtInfo.post_data < size(tempData,1)); end
                        %Remove any snippets that don't have data at the beginning of the file
                        if ~isempty(threshold_crossings); threshold_crossings = threshold_crossings(threshold_crossings - filtInfo.pre_data > 0); end
                        if ~isempty(threshold_crossings)
                            curr_spike_snippets = tempData(threshold_crossings + (-filtInfo.pre_data:filtInfo.post_data), ch);
                            curr_spike_snippets = reshape(curr_spike_snippets,[],filtInfo.snippet_length);
                            if filtInfo.align_spikes
                                peak_times = find_waveform_peaks(curr_spike_snippets', filtInfo.pre_data, filtInfo.peak_window, filtInfo.interp_factor);
                                adj_threshold_crossings = threshold_crossings + floor(peak_times)';
                                throwout_spikes = ~( adj_threshold_crossings>(filtInfo.pre_data+filtInfo.peak_offset) & (adj_threshold_crossings+(filtInfo.post_data-filtInfo.peak_offset))<size(tempData,1) );
                                threshold_crossings = threshold_crossings(~throwout_spikes);
                                adj_threshold_crossings = adj_threshold_crossings(~throwout_spikes);
                                peak_times_fraction = mod(peak_times,1);
                                curr_spike_snippets = tempData(adj_threshold_crossings + (-(filtInfo.pre_data+filtInfo.peak_offset):(filtInfo.post_data-filtInfo.peak_offset+1)), ch);
                                curr_spike_snippets = reshape(curr_spike_snippets,[],filtInfo.snippet_length+1);
                                for n = 1:size(curr_spike_snippets,1)
                                    curr_spike_snippets(n,1:(end-1)) = interp1(1:(filtInfo.snippet_length+1), curr_spike_snippets(n,:), (1:filtInfo.snippet_length)+ peak_times_fraction(n), 'spline');
                                end
                                curr_spike_snippets = curr_spike_snippets(:,1:(end-1));
                            end
                            %Note, spike times are saved based on the original threshold crossing, not adjusted for the spike alignment
                            spike_times{ch}     = [spike_times{ch}; dataBlocks.block_start_samp(tr)+uint32(threshold_crossings)-1];
                            spike_snippets{ch}  = [spike_snippets{ch}, curr_spike_snippets'];

                        end
                    end
                end
            end
        end



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

    nexFileData = fix_empty_units(nexFileData);

    %     save('A_20231005_ThresholdedData.mat', 'nexFileData', '-v7.3');  %For debugging
    writeNex5File(nexFileData, [dataPaths.save_path, envInfo.output_names{iFile} '5']);
    %     writeNexFile(nexFileData, [dataPaths.save_path, envInfo.output_names{iFile}] );
end