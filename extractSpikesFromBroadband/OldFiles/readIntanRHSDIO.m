function [intCodes,eventTimes,clockrate] = readIntanRHSDIO(dataPaths)
%read and sort intan recording file
addpath(genpath('C:\IntanDataReading\'))
[~, frequency_parameters] = read_Intan_RHS2000_file_with_path(dataPaths.input_file_path);%read_Intan_RHS2000_file_with_path(data_paths.input_file_path)
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
for k = 1:6
DIO_fileinfo = dir([dataPaths.input_file_path 'board-DIGITAL-IN-' num2str(k+8,'%02d') '.dat']);
DIO_num_samples = DIO_fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen([dataPaths.input_file_path 'board-DIGITAL-IN-' num2str(k+8,'%02d') '.dat'], 'r');
dio{k} = fread(fid, [1 DIO_num_samples], 'uint16'); %not sure about the data order here
fclose(fid);
end

clockrate = frequency_parameters.amplifier_sample_rate;
%*****Compile DIO files into discrete components
%Timestamps: Times at which events occurred 
%Events: Changes in individual bits
%Bits: Which individual bits changed
allTimestamps = [];
allEvents = [];
allBits = [];
for k=1:6
    allTimestamps = [allTimestamps; timestamps];
    allEvents = [allEvents; dio{k}'];%if digital_word = [6 num_samples]
    allBits = [allBits; (k*ones(size(dio{k})))'];
end
[allTimestamps, sort_index] = sort(allTimestamps);
allEvents = allEvents(sort_index);
allBits = allBits(sort_index);


% %subtract timestamp at creation to try to see if it fixes alignment
% allTimestamps = allTimestamps - dio{1}.first_timestamp;%timestamp_at_creation
%*****Construct matrix of discrete event codes
codeNum = 1;
eventCode = [0 0 0 0 0];
firstSix = 0;
codeInd = 1;
for n = 1 : length(allEvents)
    if allBits(n) ~= 6
        eventCode(allBits(n)) = allEvents(n);
    else
        if allEvents(n) == 1 && firstSix == 0
            firstSix = 1;
        elseif allEvents(n) == 0 && firstSix == 1
            firstSix = 0;
            if sum(eventCode) ~= 0
                codes{codeInd} = eventCode;
                intCodes(codeInd,1) = sum(eventCode .* 2.^(0:4));
                eventTimes(codeInd,1) = allTimestamps(n);
                codeInd = codeInd + 1;
            end
        end
    end
end