function calculate_Medians(envInfo, data_paths, strobeInfo, filtInfo) 

sigQual_name = [envInfo.rec_file_name, '_SigQual.mat'];
if exist([data_paths.median_path sigQual_name],'file')
    load([data_paths.median_path sigQual_name],'SigQuality')
else
    SigQuality(1).channels = [];
end

if filtInfo.overwrite_median 
    calc_median_flag = true;
else
    calc_median_flag = false;
    for ar = 1:length(envInfo.channels_to_read_by_array)
        if any(arrayfun(@(x) isequal(x.channels,envInfo.channels_to_read_by_array{ar}), SigQuality))
            calc_median_flag = true;
        end
    end

end
if calc_median_flag
 if ~isempty(strobeInfo) && ~isfield(strobeInfo, 'trial_end_samp')
        eventInfo   = openNEV( [data_paths.input_file_path , envInfo.nev_file_name], 'noread');
        %% Event strobe extracting part.
        Events = eventInfo.Data.SerialDigitalIO.UnparsedData;
        EventTimeSamples    = eventInfo.Data.SerialDigitalIO.TimeStamp;
        strobeInfo = find_trial_samples(strobeInfo, Events, EventTimeSamples);
 end
 
 %% Current function is to calculate median values.
    data_strut = readTrodesExtractedDataFile( [data_paths.input_file_path , envInfo.rec_file_name, '.raw\' envInfo.rec_file_name '.raw_nt1ch1.dat']); % %neural data
    numDataPoints = size(data_strut.fields.data,1);
    if ~isempty(strobeInfo) && length(strobeInfo.TrialIDs) > filtInfo.num_trials_for_median 
        TestTrialIDs = strobeInfo.TrialIDs(round(linspace(1, length(strobeInfo.TrialIDs), filtInfo.num_trials_for_median)));
    else
        strobeInfo.TrialIDs = 1:filtInfo.num_trials_for_median;
        TestTrialIDs = strobeInfo.TrialIDs;
        numTrialSamples = round((filtInfo.median_window/1000)*data_strut.clockrate);
        numBetweenTrialSamples = round(numDataPoints/(filtInfo.num_trials_for_median+2));
        strobeInfo.trial_start_samp = numBetweenTrialSamples + (0:(filtInfo.num_trials_for_median-1))*(numBetweenTrialSamples+numTrialSamples);
        strobeInfo.trial_end_samp = strobeInfo.trial_start_samp + numTrialSamples;
    end
    
for ar = 1:length(envInfo.channels_to_read_by_array)
    curr_index = find(arrayfun(@(x) isequal(x.channels,envInfo.channels_to_read_by_array{ar}), SigQuality));
    if  isempty(curr_index)
        if isempty(SigQuality(1).channels)
            curr_index = 1;
        else
            curr_index = length(SigQuality)+1;
        end
    elseif ~filtInfo.overwrite_median
        curr_index = [];
    end
       
    if ~isempty(curr_index)
    Data = cell(length(TestTrialIDs),1);
            
           
            SigQuality(curr_index).channels = sort(envInfo.channels_to_read_by_array{ar},'ascend');
            SigQuality(curr_index).file_channels = sort(envInfo.file_channels_to_read_by_array{ar},'ascend');
            
            for ch = 1:length(SigQuality(curr_index).channels)
                tempData = readTrodesExtractedDataFile( [data_paths.input_file_path , envInfo.rec_file_name, '.raw\' envInfo.rec_file_name '.raw_nt' num2str(SigQuality(curr_index).channels(ch)) 'ch1.dat']); 
                    
            for tr = 1:length(TestTrialIDs)
                if strobeInfo.trial_end_samp(strobeInfo.TrialIDs==TestTrialIDs(tr))<numDataPoints
                    
                    curr_time_points =  [strobeInfo.trial_start_samp(strobeInfo.TrialIDs==TestTrialIDs(tr)), strobeInfo.trial_end_samp(strobeInfo.TrialIDs==TestTrialIDs(tr))];
                      Data{tr}(:,ch) = double(tempData.fields.data(curr_time_points(1):curr_time_points(2)));
                   
                    
                    
                end
            end
            clear tempData
            end
            for tr = 1:length(TestTrialIDs)
             Data{tr} = filter(filtInfo.b, filtInfo.a, Data{tr});
            end
            SigQuality(curr_index).TestMedians = cellfun(@(x) median(abs(x)), Data, 'Uni', false);
            SigQuality(curr_index).TestMedians = cat(1,SigQuality(curr_index).TestMedians{:});
            blank_time_periods = all(SigQuality(curr_index).TestMedians<1,2);
            SigQuality(curr_index).TestMedians(SigQuality(curr_index).TestMedians<1) = NaN;
            SigQuality(curr_index).ChMedians = nanmedian(SigQuality(curr_index).TestMedians,1);
            %Check if any channels have no signal, all values equal 0 (<1 to avoid rounding errors)
            SigQuality(curr_index).BlankChannels = SigQuality(curr_index).ChMedians<1;
            SigQuality(curr_index).ChMedians(SigQuality(curr_index).BlankChannels) = 0;  %Set channels with no signal to exactly 0
            
            Data = Data(~blank_time_periods);
            tempData = cat(1,Data{:});
            tempData = tempData(:,~SigQuality(curr_index).BlankChannels);
            SigQuality(curr_index).covMat = cov(tempData);
            [SigQuality(curr_index).E, SigQuality(curr_index).D] = svd( SigQuality(curr_index).covMat );
            SigQuality(curr_index).D = diag(SigQuality(curr_index).D);
            small_val = 1e-6;
            WhtMat = SigQuality(curr_index).E * diag(1./(SigQuality(curr_index).D + small_val).^.5) * SigQuality(curr_index).E';
            WhtMat = WhtMat*(eye(size(SigQuality(curr_index).covMat)).*sqrt(SigQuality(curr_index).covMat));
            SigQuality(curr_index).WhtMat = zeros(length(SigQuality(curr_index).file_channels));
            SigQuality(curr_index).WhtMat(~SigQuality(curr_index).BlankChannels, ~SigQuality(curr_index).BlankChannels) = WhtMat;
            WhtData = cellfun(@(x) x*SigQuality(curr_index).WhtMat, Data, 'Uni', false);
            SigQuality(curr_index).WhtTestMedians = cellfun(@(x) median(abs(x)), WhtData, 'Uni', false);
            SigQuality(curr_index).WhtTestMedians = cat(1,SigQuality(curr_index).WhtTestMedians{:});
            SigQuality(curr_index).WhtChMedians = nanmedian(SigQuality(curr_index).WhtTestMedians,1);
            
            
            clear  Data WhtData tempData WhtMat
    end
end
save([data_paths.median_path sigQual_name],'SigQuality')

end

