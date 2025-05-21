function [intCodes,eventTimes,clockrate] = readIntanRHDDIO(dataPaths)
%read and sort intan recording file
% addpath(genpath('C:\IntanDataReading\'))
[~, frequency_parameters] = read_Intan_RHD2000_file(dataPaths.input_file_path);%read_Intan_RHD2000_file(data_paths.input_file_path)
clockrate = frequency_parameters.amplifier_sample_rate;
%% read DIO- board-DIGITAL-IN followed by the channel number. 
%timestamps
fileinfo = dir([dataPaths.input_file_path 'time.dat']);
t_num_samples = fileinfo.bytes/4; % int32 = 4 bytes
t_fid = fopen([dataPaths.input_file_path 'time.dat'], 'r');
timestamps = fread(t_fid, t_num_samples, 'int32');%timestamps
fclose(t_fid);
% t = t / frequency_parameters.amplifier_sample_rate; % sample rate from header file

%DIO
% "One File Per Signal Type"
% DIO_fileinfo = dir('R:\SOM RSCH\RouseLab\DataFiles\Recorded_Data\Monkey\monk_A\COTHold2022\IntanRHS\A_COTHold_20230529_230529_131418\digitalin.dat');
% num_board_channels = 6;
% DIO_num_samples = DIO_fileinfo.bytes/2; % uint16 = 2 bytes
% fid = fopen('R:\SOM RSCH\RouseLab\DataFiles\Recorded_Data\Monkey\monk_A\COTHold2022\IntanRHS\A_COTHold_20230529_230529_131418\digitalin.dat', 'r');
% digital_word = fread(fid, [1 DIO_num_samples], 'uint16'); %not sure about the data order here
% fclose(fid);
% % decint = digital_word-2^8;

%  “One File Per Channel”
% for k = 1:6
% DIO_fileinfo = dir([dataPaths.input_file_path 'board-DIGITAL-IN-' num2str(k+8,'%02d') '.dat']);
% DIO_num_samples = DIO_fileinfo.bytes/2; % uint16 = 2 bytes
% fid = fopen([dataPaths.input_file_path 'board-DIGITAL-IN-' num2str(k+8,'%02d') '.dat'], 'r');
% dio{k} = fread(fid, [1 DIO_num_samples], 'uint16'); %not sure about the data order here
% end
DIO_fileinfo = dir([dataPaths.input_file_path 'board-DIGITAL-IN-01.dat']); %could add 9-13 to debug with DIO's for event codes
DIO_num_samples = DIO_fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen([dataPaths.input_file_path 'board-DIGITAL-IN-01.dat'], 'r');
dio = fread(fid, [1 DIO_num_samples], 'uint16'); %just taking in strobes DIO
fclose(fid);
eventSamples = find(diff(dio)==1)'; %Find when stobe bit turned from 0 to 1



%*****combine Intan Events received through UDP
txt2read = dir(fullfile(dataPaths.input_file_path, '*.txt')); %Want to make sure file name matches if more than one text file exist
tmp = readtable([txt2read.folder '\' txt2read.name]);
eventCodes = tmp.Var1;  %Integer event codes
textEventSamples = tmp.Var2;  %Sample of event codes in text file
textEventSamples = textEventSamples(eventCodes>0); %Ignore the 0 event codes
eventCodes= eventCodes(eventCodes>0);

time_diffs = eventSamples-textEventSamples(1:length(eventSamples));

time_diffs = eventSamples(1:length(textEventSamples))-textEventSamples;
% figure; plot(time_diffs)

%Create temp variables of event codes and text and electronic samples

tmp_textEventSamples = textEventSamples;
tmp_eventCodes = eventCodes;
tmp_eventSamples = eventSamples;
missed_UDP_inDIO = [];
missed_DIO_inUDP = [];
missed_DIO_eventInUDP = [];
missed = 1;
counter = 0;
while  ~isempty(missed)% length(tmp_textEventSamples)>length(eventSamples)  %isempty(missed) ||
    min_length = min([length(tmp_eventSamples),length(tmp_textEventSamples)]);
    time_diffs = tmp_eventSamples(1:min_length)-tmp_textEventSamples(1:min_length);
    missed = find(abs(time_diffs)>3000,1,'first'); %Find first time when event code is off by too big of an amount
    
    if time_diffs(missed)>0 %If time difference is positive, then there was a DIO strobe but no UDP value
        missed_DIO_inUDP = [missed_DIO_inUDP tmp_textEventSamples(missed)];
        missed_DIO_eventInUDP = [missed_DIO_eventInUDP  tmp_eventCodes(missed)];
        tmp_eventSamples = ([tmp_eventSamples(1:(missed-1)); tmp_textEventSamples(missed)-1000 ;tmp_eventSamples(missed:end)]);
        counter = counter + 1;
        if counter>1000
            warning('Stuck in loop during readIntanRHDDIO.');
        end
    else  %If time difference is negative, then there was a UDP value but no DIO strobe
        missed_UDP_inDIO = [missed_UDP_inDIO tmp_eventSamples(missed)];
        tmp_textEventSamples = [tmp_textEventSamples(1:(missed-1)); tmp_eventSamples(missed)+1000 ;tmp_textEventSamples(missed:end)];
        tmp_eventCodes = [tmp_eventCodes(1:(missed-1));0;tmp_eventCodes(missed:end)];
        counter = counter + 1;
        if counter>1000
            warning('Stuck in loop during readIntanRHDDIO.');
        end
    end
    min_length = min([length(tmp_eventSamples),length(tmp_textEventSamples)]);
    time_diffs = tmp_eventSamples(1:min_length)-tmp_textEventSamples(1:min_length);
    missed = find(abs(time_diffs)>3000,1,'first');
end
tmp_eventSamples = sort(tmp_eventSamples);
tmp_textEventSamples = sort(tmp_textEventSamples);
min_length = min([length(tmp_eventSamples),length(tmp_textEventSamples)]);
time_diffs = tmp_eventSamples(1:min_length)-tmp_textEventSamples(1:min_length);

intCodes = tmp_eventCodes;
eventTimes = tmp_eventSamples;

% %*****Compile DIO files into discrete components
% %Timestamps: Times at which events occurred 
% %Events: Changes in individual bits
% %Bits: Which individual bits changed
% allTimestamps = [];
% allEvents = [];
% allBits = [];
% for k=1:6
%     allTimestamps = [allTimestamps; timestamps];
%     allEvents = [allEvents; dio{k}'];%if digital_word = [6 num_samples]
%     allBits = [allBits; (k*ones(size(dio{k})))'];
% end
% [allTimestamps, sort_index] = sort(allTimestamps);
% allEvents = allEvents(sort_index);
% allBits = allBits(sort_index);
% 
% 
% % %subtract timestamp at creation to try to see if it fixes alignment
% % allTimestamps = allTimestamps - dio{1}.first_timestamp;%timestamp_at_creation
% %*****Construct matrix of discrete event codes
% codeNum = 1;
% eventCode = [0 0 0 0 0];
% firstSix = 0;
% codeInd = 1;
% for n = 1 : length(allEvents)
%     if allBits(n) ~= 6
%         eventCode(allBits(n)) = allEvents(n);
%     else
%         if allEvents(n) == 1 && firstSix == 0
%             firstSix = 1;
%         elseif allEvents(n) == 0 && firstSix == 1
%             firstSix = 0;
%             if sum(eventCode) ~= 0
%                 codes{codeInd} = eventCode;
%                 intCodes(codeInd,1) = sum(eventCode .* 2.^(0:4));
%                 eventTimes(codeInd,1) = allTimestamps(n);
%                 codeInd = codeInd + 1;
%             end
%         end
%     end
% end

%% debug
% my_dio = cat(1,dio{:});
% my_codes = (2.^(0:4))*my_dio(1:5,:);
% my_eventTimes = [false,(my_dio(6,2:end)==1 & my_dio(6,1:(end-1))==0)];
% my_eventCodes = my_codes(my_eventTimes);
% my_eventTimes = find(my_eventTimes);
% my_eventTimes = my_eventTimes(my_eventCodes>0);
% my_eventCodes = my_eventCodes(my_eventCodes>0);

