function calculate_MediansRHD(envInfo, dataPaths, dataBlocks, filtInfo) 
% Adam Rouse, 8/13/2024, added port for Intan system to make file handling clearer and more reliable
% This function calculates the median values for electrophysiology data collected with the Intan system, 
% applies filtering, and computes whitening matrices for signal processing.

% Define the file name for storing signal quality data
sigQual_name = [envInfo.rec_file_name, '_SigQual.mat'];

% Check if the signal quality file already exists, if so, load it
if exist([dataPaths.median_path sigQual_name],'file')
    load([dataPaths.median_path sigQual_name],'SigQuality')
else
    % If not, initialize the SigQuality structure
    SigQuality(1).channels = [];
end

% Check if we should overwrite the existing median data or not
if filtInfo.overwrite_median 
    calc_median_flag = true;
else
    calc_median_flag = false;
    % Loop over each array to determine if any channel sets need re-calculation of medians
    for ar = 1:length(envInfo.channels_to_read_by_array)
        if any(arrayfun(@(x) isequal(x.channels, envInfo.channels_to_read_by_array{ar}), SigQuality))
            calc_median_flag = true;
        end
    end
end

% If re-calculation is necessary, proceed to calculate medians
if calc_median_flag

    %% Reading the data from the input file
    data_strut = readRHDExtractedDataFile([dataPaths.input_file_path 'amp-A-000.dat']);
    numDataPoints = length(data_strut.fields.data);  % Total number of data points in the file

    % Gather all .dat file names in the input path
    file_names = dir([dataPaths.input_file_path,'*.dat']);

    % If dataBlocks are not defined, initialize them using the number of trials and median window
    if isempty(dataBlocks)
        dataBlocks.BlockIDs = 1:filtInfo.num_trials_for_median;
        TestBlockIDs = dataBlocks.BlockIDs;
        numTrialSamples = round((filtInfo.median_window / 1000) * data_strut.clockrate);  % Calculate the number of samples per trial
        numBetweenTrialSamples = round((numDataPoints - numTrialSamples) / (filtInfo.num_trials_for_median + 1));  % Samples between trials
        dataBlocks.block_start_samp = numBetweenTrialSamples + (0:(filtInfo.num_trials_for_median - 1)) * numBetweenTrialSamples;
        dataBlocks.block_end_samp = dataBlocks.block_start_samp + numTrialSamples;
    end

    % Loop through each array of channels for median calculation
    for ar = 1:length(envInfo.channels_to_read_by_array)
        % Display channel number for debugging purposes
        disp(['Array ' num2str(ar) '/' num2str(length(envInfo.channels_to_read_by_array))])

        curr_index = find(arrayfun(@(x) isequal(x.channels, envInfo.channels_to_read_by_array{ar}), SigQuality));
        if isempty(curr_index)
            % If no previous signal quality data exists, create a new entry
            if isempty(SigQuality(1).channels)
                curr_index = 1;
            else
                curr_index = length(SigQuality) + 1;
            end
        elseif ~filtInfo.overwrite_median
            curr_index = [];
        end

        % Proceed if median data needs to be calculated
        if ~isempty(curr_index)
            Data = cell(length(TestBlockIDs), 1);

            % Sort and assign channel information to the SigQuality structure
            SigQuality(curr_index).channels = sort(envInfo.channels_to_read_by_array{ar}, 'ascend');
            SigQuality(curr_index).file_channels = sort(envInfo.channels_to_read_by_array{ar}, 'ascend');

            % Loop through each channel and read data from the corresponding file
            for ch = 1:length(SigQuality(curr_index).channels)
                % Display channel number for debugging purposes
                disp(['Channel ' num2str(ch) '/' num2str(length(SigQuality(curr_index).channels))])

                % Determine the corresponding data port and offset based on channel number
                if SigQuality(curr_index).channels(ch) <= 128
                    SigQuality(curr_index).port(ch) = 'A';
                    ch_offset = 0;
                elseif SigQuality(curr_index).channels(ch) > 128 && SigQuality(curr_index).channels(ch) <= 256
                    SigQuality(curr_index).port(ch) = 'B';
                    ch_offset = 128;
                elseif SigQuality(curr_index).channels(ch) <= 384
                    SigQuality(curr_index).port(ch) = 'C';
                    ch_offset = 256;
                end

                % Locate and read the corresponding .dat file for the channel
                curr_file_index = find(arrayfun(@(x) strcmpi(x.name, ['amp-' SigQuality(curr_index).port(ch) '-' num2str(SigQuality(curr_index).channels(ch) - ch_offset - 1,'%03i') '.dat']), file_names));
                if length(curr_file_index) == 1
                    tempData = readRHDExtractedDataFile([file_names(curr_file_index).folder, '\', file_names(curr_file_index).name]);

                    % Extract the data for each test block (trial)
                    for tr = 1:length(TestBlockIDs)
                        if dataBlocks.block_end_samp(dataBlocks.BlockIDs == TestBlockIDs(tr)) < numDataPoints
                            curr_time_points = [dataBlocks.block_start_samp(dataBlocks.BlockIDs == TestBlockIDs(tr)), dataBlocks.block_end_samp(dataBlocks.BlockIDs == TestBlockIDs(tr))];
                            Data{tr}(:, ch) = double(tempData.fields.data(curr_time_points(1):curr_time_points(2)));
                        end
                    end
                    clear tempData
                end
            end

            % Apply filtering to the data
            for tr = 1:length(TestBlockIDs)
                Data{tr} = filter(filtInfo.b, filtInfo.a, Data{tr});
            end

            % Compute the median values for each test block
            SigQuality(curr_index).TestMedians = cellfun(@(x) median(abs(x)), Data, 'Uni', false);
            SigQuality(curr_index).TestMedians = cat(1, SigQuality(curr_index).TestMedians{:});
            
            % Handle channels with no signal (blank channels)
            blank_time_periods = all(SigQuality(curr_index).TestMedians < 1, 2);
            SigQuality(curr_index).TestMedians(SigQuality(curr_index).TestMedians < 1) = NaN;
            SigQuality(curr_index).ChMedians = nanmedian(SigQuality(curr_index).TestMedians, 1);
            SigQuality(curr_index).BlankChannels = SigQuality(curr_index).ChMedians < 1;
            SigQuality(curr_index).ChMedians(SigQuality(curr_index).BlankChannels) = 0;  % Set channels with no signal to 0

            % Whitening process: calculate the covariance matrix and whitening transformation
            Data = Data(~blank_time_periods);
            tempData = cat(1, Data{:});
            tempData = tempData(:, ~SigQuality(curr_index).BlankChannels);
            SigQuality(curr_index).covMat = cov(tempData);
            [SigQuality(curr_index).E, SigQuality(curr_index).D] = svd(SigQuality(curr_index).covMat);
            SigQuality(curr_index).D = diag(SigQuality(curr_index).D);

            % Create whitening matrix using the covariance matrix and its singular values
            small_val = 1e-6;
            WhtMat = SigQuality(curr_index).E * diag(1 ./ (SigQuality(curr_index).D + small_val).^.5) * SigQuality(curr_index).E';
            WhtMat = WhtMat * (eye(size(SigQuality(curr_index).covMat)) .* sqrt(SigQuality(curr_index).covMat));
            SigQuality(curr_index).WhtMat = zeros(length(SigQuality(curr_index).file_channels));
            SigQuality(curr_index).WhtMat(~SigQuality(curr_index).BlankChannels, ~SigQuality(curr_index).BlankChannels) = WhtMat;

            % Apply the whitening matrix to the data
            WhtData = cellfun(@(x) x * SigQuality(curr_index).WhtMat, Data, 'Uni', false);

            % Compute the medians after whitening
            SigQuality(curr_index).WhtTestMedians = cellfun(@(x) median(abs(x)), WhtData, 'Uni', false);
            SigQuality(curr_index).WhtTestMedians = cat(1, SigQuality(curr_index).WhtTestMedians{:});
            SigQuality(curr_index).WhtChMedians = nanmedian(SigQuality(curr_index).WhtTestMedians, 1);

            % Clear temporary variables to save memory
            clear Data WhtData tempData WhtMat
        end
    end

    % Save the updated SigQuality data to the file
    save([dataPaths.median_path sigQual_name], 'SigQuality')

end
