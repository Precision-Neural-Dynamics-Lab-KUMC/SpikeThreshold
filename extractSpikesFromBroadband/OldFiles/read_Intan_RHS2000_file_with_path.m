function [amplifier_channels, frequency_parameters] = read_Intan_RHS2000_file_with_path(dataPaths)

% read_Intan_RHS2000_file
%
% Version 3.0, 8 February 2021
%
% Reads Intan Technologies RHS data file generated by Intan Stimulation /
% Recording Controller.  Data are parsed and placed into variables that
% appear in the base MATLAB workspace.  Therefore, it is recommended to
% execute a 'clear' command before running this program to clear all other
% variables from the base workspace.
%
% Example:
% >> clear
% >> read_Intan_RHS2000_file
% >> whos
% >> amplifier_channels(1)
% >> plot(t, amplifier_data(1,:))

% [file, path, filterindex] = ...
%     uigetfile('*.rhs', 'Select an RHS2000 Data File', 'MultiSelect', 'off');
% curr_monkey
% if ~exist('dataPaths', 'var')
%     recorded_data_location = ['\\kumc.edu\data\research\SOM RSCH\RouseLab\DataFiles\Recorded_Data\Monkey\monk_' curr_monkey, '\', curr_task];
% end
fileList = dir(fullfile(dataPaths,'*.rhs'));
if (fileList.name == 0)
    return;
end

tic;
filename = [fileList.folder '\' fileList.name];
fid = fopen(filename, 'r');

s = dir(filename);
filesize = s.bytes;

% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHS2000 data file.
magic_number = fread(fid, 1, 'uint32');
if magic_number ~= hex2dec('d69127ac')
    error('Unrecognized file type.');
end

% Read version number.
data_file_main_version_number = fread(fid, 1, 'int16');
data_file_secondary_version_number = fread(fid, 1, 'int16');

% fprintf(1, '\n');
% fprintf(1, 'Reading Intan Technologies RHS2000 Data File, Version %d.%d\n', ...
%     data_file_main_version_number, data_file_secondary_version_number);
% fprintf(1, '\n');
% 
num_samples_per_data_block = 128;

% Read information of sampling rate and amplifier frequency settings.
sample_rate = fread(fid, 1, 'single');
dsp_enabled = fread(fid, 1, 'int16');
actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
actual_lower_bandwidth = fread(fid, 1, 'single');
actual_lower_settle_bandwidth = fread(fid, 1, 'single');
actual_upper_bandwidth = fread(fid, 1, 'single');

desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
desired_lower_bandwidth = fread(fid, 1, 'single');
desired_lower_settle_bandwidth = fread(fid, 1, 'single');
desired_upper_bandwidth = fread(fid, 1, 'single');

% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
notch_filter_mode = fread(fid, 1, 'int16');
notch_filter_frequency = 0;
if (notch_filter_mode == 1)
    notch_filter_frequency = 50;
elseif (notch_filter_mode == 2)
    notch_filter_frequency = 60;
end

desired_impedance_test_frequency = fread(fid, 1, 'single');
actual_impedance_test_frequency = fread(fid, 1, 'single');

amp_settle_mode = fread(fid, 1, 'int16');
charge_recovery_mode = fread(fid, 1, 'int16');

stim_step_size = fread(fid, 1, 'single');
charge_recovery_current_limit = fread(fid, 1, 'single');
charge_recovery_target_voltage = fread(fid, 1, 'single');

% Place notes in data strucure
notes = struct( ...
    'note1', fread_QString(fid), ...
    'note2', fread_QString(fid), ...
    'note3', fread_QString(fid) );
    
% See if dc amplifier data was saved
dc_amp_data_saved = fread(fid, 1, 'int16');

% Load board mode.
board_mode = fread(fid, 1, 'int16');

reference_channel = fread_QString(fid);

% Place frequency-related information in data structure.
frequency_parameters = struct( ...
    'amplifier_sample_rate', sample_rate, ...
    'board_adc_sample_rate', sample_rate, ...
    'board_dig_in_sample_rate', sample_rate, ...
    'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
    'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
    'dsp_enabled', dsp_enabled, ...
    'desired_lower_bandwidth', desired_lower_bandwidth, ...
    'desired_lower_settle_bandwidth', desired_lower_settle_bandwidth, ...
    'actual_lower_bandwidth', actual_lower_bandwidth, ...
    'actual_lower_settle_bandwidth', actual_lower_settle_bandwidth, ...
    'desired_upper_bandwidth', desired_upper_bandwidth, ...
    'actual_upper_bandwidth', actual_upper_bandwidth, ...
    'notch_filter_frequency', notch_filter_frequency, ...
    'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
    'actual_impedance_test_frequency', actual_impedance_test_frequency );

stim_parameters = struct( ...
    'stim_step_size', stim_step_size, ...
    'charge_recovery_current_limit', charge_recovery_current_limit, ...
    'charge_recovery_target_voltage', charge_recovery_target_voltage, ...
    'amp_settle_mode', amp_settle_mode, ...
    'charge_recovery_mode', charge_recovery_mode );

% Define data structure for spike trigger settings.
spike_trigger_struct = struct( ...
    'voltage_trigger_mode', {}, ...
    'voltage_threshold', {}, ...
    'digital_trigger_channel', {}, ...
    'digital_edge_polarity', {} );

new_trigger_channel = struct(spike_trigger_struct);
spike_triggers = struct(spike_trigger_struct);

% Define data structure for data channels.
channel_struct = struct( ...
    'native_channel_name', {}, ...
    'custom_channel_name', {}, ...
    'native_order', {}, ...
    'custom_order', {}, ...
    'board_stream', {}, ...
    'chip_channel', {}, ...
    'port_name', {}, ...
    'port_prefix', {}, ...
    'port_number', {}, ...
    'electrode_impedance_magnitude', {}, ...
    'electrode_impedance_phase', {} );

new_channel = struct(channel_struct);

% Create structure arrays for each type of data channel.
amplifier_channels = struct(channel_struct);
board_adc_channels = struct(channel_struct);
board_dac_channels = struct(channel_struct);
board_dig_in_channels = struct(channel_struct);
board_dig_out_channels = struct(channel_struct);

amplifier_index = 1;
board_adc_index = 1;
board_dac_index = 1;
board_dig_in_index = 1;
board_dig_out_index = 1;

% Read signal summary from data file header.

number_of_signal_groups = fread(fid, 1, 'int16');

for signal_group = 1:number_of_signal_groups
    signal_group_name = fread_QString(fid);
    signal_group_prefix = fread_QString(fid);
    signal_group_enabled = fread(fid, 1, 'int16');
    signal_group_num_channels = fread(fid, 1, 'int16');
    signal_group_num_amp_channels = fread(fid, 1, 'int16');

    if (signal_group_num_channels > 0 && signal_group_enabled > 0)
        new_channel(1).port_name = signal_group_name;
        new_channel(1).port_prefix = signal_group_prefix;
        new_channel(1).port_number = signal_group;
        for signal_channel = 1:signal_group_num_channels
            new_channel(1).native_channel_name = fread_QString(fid);
            new_channel(1).custom_channel_name = fread_QString(fid);
            new_channel(1).native_order = fread(fid, 1, 'int16');
            new_channel(1).custom_order = fread(fid, 1, 'int16');
            signal_type = fread(fid, 1, 'int16');
            channel_enabled = fread(fid, 1, 'int16');
            new_channel(1).chip_channel = fread(fid, 1, 'int16');
            fread(fid, 1, 'int16');  % ignore command_stream
            new_channel(1).board_stream = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
            new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
            new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
            
            if (channel_enabled)
                switch (signal_type)
                    case 0
                        amplifier_channels(amplifier_index) = new_channel;
                        spike_triggers(amplifier_index) = new_trigger_channel;
                        amplifier_index = amplifier_index + 1;
                    case 1
                        % aux inputs; not used in RHS2000 system
                    case 2
                        % supply voltage; not used in RHS2000 system
                    case 3
                        board_adc_channels(board_adc_index) = new_channel;
                        board_adc_index = board_adc_index + 1;
                    case 4
                        board_dac_channels(board_dac_index) = new_channel;
                        board_dac_index = board_dac_index + 1;
                    case 5
                        board_dig_in_channels(board_dig_in_index) = new_channel;
                        board_dig_in_index = board_dig_in_index + 1;
                    case 6
                        board_dig_out_channels(board_dig_out_index) = new_channel;
                        board_dig_out_index = board_dig_out_index + 1;
                    otherwise
                        error('Unknown channel type');
                end
            end
            
        end
    end
end

% Summarize contents of data file.
num_amplifier_channels = amplifier_index - 1;
num_board_adc_channels = board_adc_index - 1;
num_board_dac_channels = board_dac_index - 1;
num_board_dig_in_channels = board_dig_in_index - 1;
num_board_dig_out_channels = board_dig_out_index - 1;

% fprintf(1, 'Found %d amplifier channel%s.\n', ...
%     num_amplifier_channels, plural(num_amplifier_channels));
% if (dc_amp_data_saved ~= 0)
%     fprintf(1, 'Found %d DC amplifier channel%s.\n', ...
%         num_amplifier_channels, plural(num_amplifier_channels));
% end
% fprintf(1, 'Found %d board ADC channel%s.\n', ...
%     num_board_adc_channels, plural(num_board_adc_channels));
% fprintf(1, 'Found %d board DAC channel%s.\n', ...
%     num_board_dac_channels, plural(num_board_dac_channels));
% fprintf(1, 'Found %d board digital input channel%s.\n', ...
%     num_board_dig_in_channels, plural(num_board_dig_in_channels));
% fprintf(1, 'Found %d board digital output channel%s.\n', ...
%     num_board_dig_out_channels, plural(num_board_dig_out_channels));
% fprintf(1, '\n');
% 
% Determine how many samples the data file contains.

% Each data block contains num_samples_per_data_block amplifier samples.
bytes_per_block = num_samples_per_data_block * 4;  % timestamp data
if (dc_amp_data_saved ~= 0)
    bytes_per_block = bytes_per_block + num_samples_per_data_block * (2 + 2 + 2) * num_amplifier_channels;
else
    bytes_per_block = bytes_per_block + num_samples_per_data_block * (2 + 2) * num_amplifier_channels;    
end
% Board analog inputs are sampled at same rate as amplifiers
bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_adc_channels;
% Board analog outputs are sampled at same rate as amplifiers
bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_dac_channels;
% Board digital inputs are sampled at same rate as amplifiers
if (num_board_dig_in_channels > 0)
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
end
% Board digital outputs are sampled at same rate as amplifiers
if (num_board_dig_out_channels > 0)
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
end

% How many data blocks remain in this file?
data_present = 0;
bytes_remaining = filesize - ftell(fid);
if (bytes_remaining > 0)
    data_present = 1;
end

num_data_blocks = bytes_remaining / bytes_per_block;

num_amplifier_samples = num_samples_per_data_block * num_data_blocks;
num_board_adc_samples = num_samples_per_data_block * num_data_blocks;
num_board_dac_samples = num_samples_per_data_block * num_data_blocks;
num_board_dig_in_samples = num_samples_per_data_block * num_data_blocks;
num_board_dig_out_samples = num_samples_per_data_block * num_data_blocks;

record_time = num_amplifier_samples / sample_rate;

% if (data_present)
%     fprintf(1, 'File contains %0.3f seconds of data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
%         record_time, sample_rate / 1000);
%     fprintf(1, '\n');
% else
%     fprintf(1, 'Header file contains no data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
%         sample_rate / 1000);
%     fprintf(1, '\n');
% end

if (data_present)
    
    % Pre-allocate memory for data.
    fprintf(1, 'Allocating memory for data...\n');

    t = zeros(1, num_amplifier_samples);

    amplifier_data = zeros(num_amplifier_channels, num_amplifier_samples);
    if (dc_amp_data_saved ~= 0)
        dc_amplifier_data = zeros(num_amplifier_channels, num_amplifier_samples);
    end
    stim_data = zeros(num_amplifier_channels, num_amplifier_samples);
    amp_settle_data = zeros(num_amplifier_channels, num_amplifier_samples);
    charge_recovery_data = zeros(num_amplifier_channels, num_amplifier_samples);
    compliance_limit_data = zeros(num_amplifier_channels, num_amplifier_samples);
    board_adc_data = zeros(num_board_adc_channels, num_board_adc_samples);
    board_dac_data = zeros(num_board_dac_channels, num_board_dac_samples);
    board_dig_in_data = zeros(num_board_dig_in_channels, num_board_dig_in_samples);
    board_dig_in_raw = zeros(1, num_board_dig_in_samples);
    board_dig_out_data = zeros(num_board_dig_out_channels, num_board_dig_out_samples);
    board_dig_out_raw = zeros(1, num_board_dig_out_samples);

    % Read sampled data from file.
    fprintf(1, 'Reading data from file...\n');

    amplifier_index = 1;
    board_adc_index = 1;
    board_dac_index = 1;
    board_dig_in_index = 1;
    board_dig_out_index = 1;

    print_increment = 10;
    percent_done = print_increment;
    for i=1:num_data_blocks
        t(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'int32');
        if (num_amplifier_channels > 0)
            amplifier_data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
            if (dc_amp_data_saved ~= 0)
                dc_amplifier_data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
            end
            stim_data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
        end
        if (num_board_adc_channels > 0)
            board_adc_data(:, board_adc_index:(board_adc_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_board_adc_channels], 'uint16')';
        end
        if (num_board_dac_channels > 0)
            board_dac_data(:, board_dac_index:(board_dac_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_board_dac_channels], 'uint16')';
        end
        if (num_board_dig_in_channels > 0)
            board_dig_in_raw(board_dig_in_index:(board_dig_in_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
        end
        if (num_board_dig_out_channels > 0)
            board_dig_out_raw(board_dig_out_index:(board_dig_out_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
        end

        amplifier_index = amplifier_index + num_samples_per_data_block;
        board_adc_index = board_adc_index + num_samples_per_data_block;
        board_dac_index = board_dac_index + num_samples_per_data_block;
        board_dig_in_index = board_dig_in_index + num_samples_per_data_block;
        board_dig_out_index = board_dig_out_index + num_samples_per_data_block;

        fraction_done = 100 * (i / num_data_blocks);
        if (fraction_done >= percent_done)
            fprintf(1, '%d%% done...\n', percent_done);
            percent_done = percent_done + print_increment;
        end
    end

    % Make sure we have read exactly the right amount of data.
    bytes_remaining = filesize - ftell(fid);
    if (bytes_remaining ~= 0)
        %error('Error: End of file not reached.');
    end

end

% Close data file.
fclose(fid);

if (data_present)
    
    fprintf(1, 'Parsing data...\n');

    % Extract digital input channels to separate variables.
    for i=1:num_board_dig_in_channels
       mask = 2^(board_dig_in_channels(i).native_order) * ones(size(board_dig_in_raw));
       board_dig_in_data(i, :) = (bitand(board_dig_in_raw, mask) > 0);
    end
    for i=1:num_board_dig_out_channels
       mask = 2^(board_dig_out_channels(i).native_order) * ones(size(board_dig_out_raw));
       board_dig_out_data(i, :) = (bitand(board_dig_out_raw, mask) > 0);
    end

    % Scale voltage levels appropriately.
    amplifier_data = 0.195 * (amplifier_data - 32768); % units = microvolts
    if (dc_amp_data_saved ~= 0)
        dc_amplifier_data = -0.01923 * (dc_amplifier_data - 512); % units = volts
    end
    compliance_limit_data = stim_data >= 2^15;
    stim_data = stim_data - (compliance_limit_data * 2^15);
    charge_recovery_data = stim_data >= 2^14;
    stim_data = stim_data - (charge_recovery_data * 2^14);
    amp_settle_data = stim_data >= 2^13;
    stim_data = stim_data - (amp_settle_data * 2^13);
    stim_polarity = stim_data >= 2^8;
    stim_data = stim_data - (stim_polarity * 2^8);
    stim_polarity = 1 - 2 * stim_polarity; % convert (0 = pos, 1 = neg) to +/-1
    stim_data = stim_data .* stim_polarity;
    stim_data = stim_parameters.stim_step_size * stim_data / 1.0e-6; % units = microamps
    board_adc_data = 312.5e-6 * (board_adc_data - 32768); % units = volts
    board_dac_data = 312.5e-6 * (board_dac_data - 32768); % units = volts

    % Check for gaps in timestamps.
    num_gaps = sum(diff(t) ~= 1);
    if (num_gaps == 0)
        fprintf(1, 'No missing timestamps in data.\n');
    else
        fprintf(1, 'Warning: %d gaps in timestamp data found.  Time scale will not be uniform!\n', ...
            num_gaps);
    end

    % Scale time steps (units = seconds).
    t = t / sample_rate;

    % If the software notch filter was selected during the recording, apply the
    % same notch filter to amplifier data here.  But don't do this for v3.0+ 
    % files (from Intan RHX software) because RHX saves notch-filtered data.
    if (notch_filter_frequency > 0 && data_file_main_version_number < 3)
        fprintf(1, 'Applying notch filter...\n');

        print_increment = 10;
        percent_done = print_increment;
        for i=1:num_amplifier_channels
            amplifier_data(i,:) = ...
                notch_filter(amplifier_data(i,:), sample_rate, notch_filter_frequency, 10);

            fraction_done = 100 * (i / num_amplifier_channels);
            if (fraction_done >= percent_done)
                fprintf(1, '%d%% done...\n', percent_done);
                percent_done = percent_done + print_increment;
            end

        end
    end

end

% Move variables to base workspace.

move_to_base_workspace(notes);
move_to_base_workspace(frequency_parameters);
move_to_base_workspace(stim_parameters);
if (data_file_main_version_number > 1)
    move_to_base_workspace(reference_channel);
end

if (num_amplifier_channels > 0)
    move_to_base_workspace(amplifier_channels);
    if (data_present)
        move_to_base_workspace(amplifier_data);
        if (dc_amp_data_saved ~= 0)
            move_to_base_workspace(dc_amplifier_data);
        end
        move_to_base_workspace(stim_data);
        move_to_base_workspace(amp_settle_data);
        move_to_base_workspace(charge_recovery_data);
        move_to_base_workspace(compliance_limit_data);
        move_to_base_workspace(t);
    end
    move_to_base_workspace(spike_triggers);
end
if (num_board_adc_channels > 0)
    move_to_base_workspace(board_adc_channels);
    if (data_present)
        move_to_base_workspace(board_adc_data);
    end
end
if (num_board_dac_channels > 0)
    move_to_base_workspace(board_dac_channels);
    if (data_present)
        move_to_base_workspace(board_dac_data);
    end
end
if (num_board_dig_in_channels > 0)
    move_to_base_workspace(board_dig_in_channels);
    if (data_present)
        move_to_base_workspace(board_dig_in_data);
    end
end
if (num_board_dig_out_channels > 0)
    move_to_base_workspace(board_dig_out_channels);
    if (data_present)
        move_to_base_workspace(board_dig_out_data);
    end
end

% fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', toc);
% if (data_present)
%     fprintf(1, 'Extracted data are now available in the MATLAB workspace.\n');
% else
%     fprintf(1, 'Extracted waveform information is now available in the MATLAB workspace.\n');
% end
% fprintf(1, 'Type ''whos'' to see variables.\n');
% fprintf(1, '\n');

return


function a = fread_QString(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
% the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
% convert length from bytes to 16-bit Unicode words
length = length / 2;

for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end

return


function s = plural(n)

% s = plural(n)
% 
% Utility function to optionally plurailze words based on the value
% of n.

if (n == 1)
    s = '';
else
    s = 's';
end

return


function out = notch_filter(in, fSample, fNotch, Bandwidth)

% out = notch_filter(in, fSample, fNotch, Bandwidth)
%
% Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
% fSample = sample rate of data (in Hz or Samples/sec)
% fNotch = filter notch frequency (in Hz)
% Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
%   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
%   poor time-domain properties with an extended ringing response to
%   transient disturbances.
%
% Example:  If neural data was sampled at 30 kSamples/sec
% and you wish to implement a 60 Hz notch filter:
%
% out = notch_filter(in, 30000, 60, 10);

tstep = 1/fSample;
Fc = fNotch*tstep;

L = length(in);

% Calculate IIR filter parameters
d = exp(-2*pi*(Bandwidth/2)*tstep);
b = (1 + d*d)*cos(2*pi*Fc);
a0 = 1;
a1 = -b;
a2 = d*d;
a = (1 + d*d)/2;
b0 = 1;
b1 = -2*cos(2*pi*Fc);
b2 = 1;

out = zeros(size(in));
out(1) = in(1);  
out(2) = in(2);
% (If filtering a continuous data stream, change out(1) and out(2) to the
%  previous final two values of out.)

% Run filter
for i=3:L
    out(i) = (a*b2*in(i-2) + a*b1*in(i-1) + a*b0*in(i) - a2*out(i-2) - a1*out(i-1))/a0;
end

return


function move_to_base_workspace(variable)

% move_to_base_workspace(variable)
%
% Move variable from function workspace to base MATLAB workspace so
% user will have access to it after the program ends.

variable_name = inputname(1);
assignin('base', variable_name, variable);

return;

