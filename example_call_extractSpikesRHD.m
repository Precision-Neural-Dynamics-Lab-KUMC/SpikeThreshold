%Example call script for using extractSpikesFromBroadband package
%Change code with ******* to your path and file information


%*********Added tools needed: NPMK and HowToReadAndWriteNexAndNex5FilesInMatlab
% addpath(genpath([root_file_path 'added_matlab_tools\NPMK-4.5.3.0']))
% addpath('C:\data_processing\HowToReadAndWriteNexAndNex5FilesInMatlab')

%********Add extractSpikesFromBroadBand to path
% addpath(genpath('.\extractSpikesFromBroadBand'))

%*****File name to read
%%% Input from DailyProcessing
% filename = 'A_COTPinball_20220204';
envInfo.rec_file_name	= filename;


% %*****Currently set up assuming .ns5 and .nev have same file name, can also
% %be edited here with a name independent of .ns5 name
% envInfo.nev_file_name = regexprep( envInfo.ns5_file_name, 'ns\d', 'nev');  


%*********Each element of cell is a group of channels in an array that will have
%their covariance calculated and the common average subtracted away
%envInfo.channels_to_read_by_array = {1:32,33:64,65:96,97:128,129:160,161:192,193:224,225:256,257:288,289:320,321:352,353:384};%2 ped to analyze
envInfo.channels_to_read_by_array = ch_to_process;

%*******The output .nex file(s) are specified here
%.output_names is the .nex file names to be output
envInfo.output_names = {[filename '_out.nex']};%_out.nex : 
%.array_to_fileNum is which file in .output_names a given array in
%.channels_to_read_by_array should be output 
% All 1's simply outputs all specified channels to the first file name
envInfo.array_to_fileNum = ones(size(envInfo.channels_to_read_by_array));

% envInfo.channels_to_read_by_array = {33:48,49:64};
% envInfo.output_names = {'P_20170705_Iout.nex', 'P_20170705_J_out.nex'};
% envInfo.array_to_fileNum = [1 2 ];



 
% root_data_path = '\\kumc.edu\data\research\SOM RSCH\RouseLab\DataFiles\';
dataPaths.input_file_path = spikeThreshPaths.input_file_path;
% dataPaths.kinarm_filename = matProcessingFiles.kinarm_filename;
% dataPaths.input_file_path	= [root_data_path filename '.rec\'];
% % dataPaths.input_file_path	= 'C:\DataFiles\data_raw\monk_p\20160504_COT_precision\P_CerebusData\';
dataPaths.median_path		= [root_data_path 'Processed_Data\' curr_task '\monk_' curr_monkey '\SpikeQuality\'];
dataPaths.save_path			= [root_data_path 'Processed_Data\' curr_task '\monk_' curr_monkey '\BBtoSpikes\'];


%%% Input from DailyProcessing
% root_data_path = '\\kumc.edu\data\research\SOM RSCH\RouseLab\DataFiles\';%'E:\Data\'
% dataPaths.input_file_path	= [root_data_path 'Recorded_Data\Monkey\monk_A\COTPerturb20210713\SpikeGadgets\' filename '.rec\'];
% dataPaths.median_path		= [root_data_path 'Processed_Data\COTPerturb20210713\monk_A\SpikeQuality\'];
% dataPaths.save_path			= [root_data_path 'Processed_Data\COTPerturb20210713\monk_A\BBtoSpikes\'];


filtInfo = defaultFiltInfo();

% strobeInfo.trial_start_strb    = 'TrialID';  % Special case of TrialID uses all event codes <3000, otherwise works with a single event code for start of trial
% strobeInfo.trial_end_strb      = 6013;  % Right now, the trial_end_code must be immediately before the next trial_start_code, TODO - add more flexibility to time between end and next start code
% strobeInfo.spike_end_strb      = 6013;  % If I put the particular event strobe, then I can get part of the data
% strobeInfo.spike_end_offset    = 0;   % (value: 0 to minus integer) -integer is number of samples before spike_end_strb, putting 0 here will go all the until spike_end_strb.

% strobeInfo.trial_start_strb    = 'TrialID';  % Special case of TrialID uses all event codes <3000, otherwise works with a single event code for start of trial
% strobeInfo.trial_end_strb      = 6253;  % Right now, the trial_end_code must be immediately before the next trial_start_code, TODO - add more flexibility to time between end and next start code
% strobeInfo.spike_end_strb      = 6120;  % If I put the particular event strobe, then I can get part of the data
% strobeInfo.spike_end_offset    = -10;   % (value: 0 to minus integer) -integer is number of samples before spike_end_strb, putting 0 here will go all the until spike_end_strb.

strobeInfo = [];
extractSpikesRHD(dataPaths, envInfo, strobeInfo, filtInfo)


