function [intCodes,eventTimes,clockrate] =  readspikegadgetsDIO(data_DIO_filename)

%*****Read extracted DIO files into matlab
dioFilename = [data_DIO_filename '.dio_Controller_Din'];
for k = 1:6
dio{k} = readTrodesExtractedDataFile([dioFilename num2str(k) '.dat']);
end
clockrate = dio{1}.clockrate;
%*****Compile DIO files into discrete components
%Timestamps: Times at which events occurred 
%Events: Changes in individual bits
%Bits: Which individual bits changed
allTimestamps = [];
allEvents = [];
allBits = [];
for k = 1:6
    allTimestamps = [allTimestamps; dio{k}.fields(1).data];
    allEvents = [allEvents; dio{k}.fields(2).data];
    allBits = [allBits; k*ones(size(dio{k}.fields(2).data))];
end
[allTimestamps, sort_index] = sort(allTimestamps);
allEvents = allEvents(sort_index);
allBits = allBits(sort_index);

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
                intCodes(codeInd) = sum(eventCode .* 2.^(0:4));
                eventTimes(codeInd) = allTimestamps(n);
                codeInd = codeInd + 1;
            end
        end
    end
end
