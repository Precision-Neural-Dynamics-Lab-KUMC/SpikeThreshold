function [Events, EventTimeSamples, sampling_rate] = extractEvents(dataPaths, kinarmData)
%% extract event codes from Kinarm and Intan DIO/UDP
%Comparing KinArmevent and Intan DIO/UDP event
%Willy, Adam, Lauren, and Xavier, 9/11/24
%There are three types of event times/samples and two arrays of event codes
%Kinarm file - event time and codes
%Intan DIO - event times only
%Intan UDP - event time and codes


% Reads KinArm data and extract eventcodes and timestamp in seconds
% event_codes KinArm
% 1-'TRIAL_START'
% 2-'STAY_CENTER'
% 3-'TARGET_ON'
% 4-'HOLD_TARGET'
% 5-'REWARD_SUCCESS'
%
% 9-'LEAVE_CENTER'
% 10-'PERTURB'
% 15-'TRIAL_END'
% (99)-'TASK_RESET'


[intCodes,eventSamples,textEventSamples, sampling_rate] = readIntanRHDDIO_nocorrection(dataPaths);

KinArm_eventTimes = [];
KinArm_eventCodes = [];
KinArm_eventLabels = [];

for tr=1:length(kinarmData.c3d)
    timeStr = kinarmData.c3d(tr).TRIAL.TIME;
    timeObj = datetime(timeStr, 'Format', 'HH:mm:ss.SSSSSS');
    timeInSeconds = timeObj.Hour * 3600 + timeObj.Minute * 60 + timeObj.Second;
    KinArm_eventTimes = [KinArm_eventTimes kinarmData.c3d(tr).EVENTS.TIMES + timeInSeconds];
    KinArm_eventLabels = [KinArm_eventLabels kinarmData.c3d(tr).EVENTS.LABELS];
    %covert to 3D trial ID's
    tr_id = kinarmData.c3d(tr).TRIAL.TRIAL_NUM;
    % Compute the codes based on the transformation-tr_id = code1*16^2+code2*16+code3
    code1 = floor(tr_id / 16^2);          % Most significant part
    remaining = mod(tr_id, 16^2);         % Remaining part after extracting code1
    code2 = floor(remaining / 16);              % Middle part
    code3 = mod(remaining, 16);                 % Least significant part

    % Adjust the codes to start from the required base values
    code1 = code1 + 144;  % Adjusting base for code1
    code2 = code2 + 80;   % Adjusting base for code2
    code3 = code3 + 16;   % Adjusting base for code3

    %assigns the event codes to a pre-allocated matrix base on the event label
    tmp_eventcodes = nan(length(kinarmData.c3d(tr).EVENTS.LABELS),1);
    tmp_eventcodes(1:3) = [code1 code2 code3];%trial ID codes
    for j = 1: length(kinarmData.c3d(1).EVENT_DEFINITIONS.LABELS)
        tmp_eventcodes(strcmpi(kinarmData.c3d(tr).EVENTS.LABELS,kinarmData.c3d(1).EVENT_DEFINITIONS.LABELS{j}))=kinarmData.c3d(1).EVENT_DEFINITIONS.CODES(j);
    end
    KinArm_eventCodes = [KinArm_eventCodes; tmp_eventcodes];
end

% NaN in KinArm_eventCodes are likely 'TASK_BUTTON_X_CLICKED' which are not
% sent to intan
KinArm_eventTimes = KinArm_eventTimes(~isnan(KinArm_eventCodes));
KinArm_eventCodes = KinArm_eventCodes(~isnan(KinArm_eventCodes));
KinArm_eventTimes = KinArm_eventTimes - KinArm_eventTimes(1);
KinArm_eventTimes = KinArm_eventTimes';

%The text events have a tendency to have the same time stamp at the start
%of a trial, here we go and add a minimum difference to have each text
%event have its own unique time stamp to help with future code
min_diff_samp = round(min(diff(KinArm_eventTimes)).*30000);
while (sum(diff(textEventSamples)<min_diff_samp) ) > 0
    textEventSamples([false; diff(textEventSamples)<min_diff_samp]) = textEventSamples([diff(textEventSamples)<min_diff_samp; false]) + min_diff_samp;
end

% disp(filename)
disp(['KinArm events: ' num2str(length(KinArm_eventCodes))])
disp(['DIO events: ' num2str(length(eventSamples))])
disp(['Text events: ' num2str(length(textEventSamples))])
disp(' ')
min_length = min([length(eventSamples),length(textEventSamples)]);

% time_diffs = eventSamples(1:min_length)-textEventSamples(1:min_length);
% figure; plot(time_diffs)
% ylim([0,16e4])
% title(filename, 'Interpreter', 'none')

%Intan is usually sampled at 30 kHz (although might be 20 kHz if set incorrectly!), data originally in samples, get dio/UDP times
textEventSamplesTimes = textEventSamples./sampling_rate;

%Uses the Kinarm event codes as ground truth to line up the Intan UDP text event
%codes
correctedintcodes = nan(size(KinArm_eventCodes));
corrected_textEventTimes = nan(size(KinArm_eventCodes));
corrected_textEventSamples = nan(size(KinArm_eventCodes));
c = 1;
for k = 1: length(KinArm_eventCodes)
    if c <= length(intCodes)
        if KinArm_eventCodes(k) == intCodes(c)
            correctedintcodes(k) = intCodes(c);
            corrected_textEventTimes(k) = textEventSamplesTimes(c);
            corrected_textEventSamples(k) = textEventSamples(c);
            c = c + 1;
        end
    end
end

missedTextEvents = isnan(correctedintcodes);
disp(['Missed Text events: ' num2str(sum(missedTextEvents))])

%For missed UDP text events, make an estimate based on average time delay
%from Kinarm to Intan UDP
timedifference = corrected_textEventTimes - KinArm_eventTimes;
mediantimedifference = median(timedifference,'omitnan');
corrected_textEventTimes(missedTextEvents) = KinArm_eventTimes(missedTextEvents) + mediantimedifference;
corrected_textEventSamples(missedTextEvents) = round(corrected_textEventTimes(missedTextEvents)*sampling_rate);

% Now, compare Intan DIO and UDP,
% Find the most most common difference in samples between the events
% unclear why but DIO event samples appear to be after text event samples
sample_diffs  = eventSamples - corrected_textEventSamples';
edges = -6000:100:6000;
tmp = discretize(sample_diffs,edges,(edges(1:(end-1))+edges(2:end))./2);
best_lag = mode(tmp(:));

% figure; scatter(corrected_textEventSamples, ones(size(corrected_textEventSamples)))
% hold on; scatter(eventSamples, 1.01*ones(size(eventSamples)))
% ylim([0,2])

corrected_eventSamples = NaN(size(corrected_textEventSamples));
counter = 1;
for k = 1:length(eventSamples)
    while counter <= length(corrected_textEventSamples) && ...
            ( eventSamples(k) - corrected_textEventSamples(counter) ) > (best_lag+3000) %Keep iterating until the text event sample is close to the current DIO event sample (the best lag plus 0.1s)
        counter = counter+1;
    end
    if counter <= length(corrected_textEventSamples)
        corrected_eventSamples(counter) = eventSamples(k);
    end
    counter = counter+1;
end

missedDigEvents = isnan(corrected_eventSamples);
disp(['Missed DIO events: ' num2str(sum(missedDigEvents))])

[b,~,~,~,stats] = regress(corrected_textEventSamples(~missedTextEvents)./sampling_rate,[KinArm_eventTimes(~missedTextEvents), ones(size(KinArm_eventTimes(~missedTextEvents)))]);

corrected_textEventSamples(missedTextEvents) = round( (b(1)*KinArm_eventTimes(missedTextEvents) + b(2)) .*sampling_rate );

% This includes digital events that had a text file event but no digital and
% events only seen in the kinarm file
corrected_eventSamples(missedDigEvents) = corrected_textEventSamples(missedDigEvents)+best_lag;

% We then write back the missed text event samples from digital samples,
% should be the same if neither events were there originally but uses
% digital events if possible since more closely related than Kinarm
corrected_textEventSamples(missedTextEvents) = corrected_eventSamples(missedTextEvents)-best_lag;

textKinDiffs = corrected_textEventSamples./sampling_rate - KinArm_eventTimes;
dioKinDiffs = corrected_eventSamples./sampling_rate - KinArm_eventTimes;

figure
scatter(KinArm_eventTimes(~missedTextEvents), textKinDiffs(~missedTextEvents))
hold on
scatter(KinArm_eventTimes(missedTextEvents), textKinDiffs(missedTextEvents))
legend('Present', 'Missed', 'location', 'Southeast')
title('Kinarm times vs. UDP text times')

figure
scatter(KinArm_eventTimes(~missedDigEvents), dioKinDiffs(~missedDigEvents))
hold on
scatter(KinArm_eventTimes(missedDigEvents), dioKinDiffs(missedDigEvents))
legend('Present', 'Missed', 'location', 'Southeast')
title('Kinarm times vs. DIO times')

%perform regression only using data clearly aligned between KinArm and DIO
max_timediff = 0.001;
eventsForRegress = [diff(dioKinDiffs)<max_timediff; false] & [false; diff(dioKinDiffs)<max_timediff];

[b,~,~,~,stats] = regress(corrected_eventSamples(~missedDigEvents&eventsForRegress)./sampling_rate,[KinArm_eventTimes(~missedDigEvents&eventsForRegress), ones(size(KinArm_eventTimes(~missedDigEvents&eventsForRegress)))]);

pred_eventSamples = round( (b(1)*KinArm_eventTimes + b(2)) .*sampling_rate );

predEventDiff = corrected_eventSamples - pred_eventSamples;

%Adjust all the missed digital events to the predicted value
corrected_eventSamples(missedDigEvents) = pred_eventSamples(missedDigEvents);
max_allowed_diff = 50;
eventsToAdjust = abs(predEventDiff)>max_allowed_diff;
%Adjust outliers to predict based off Kinarm as that's likely more accurate
corrected_eventSamples(eventsToAdjust) = pred_eventSamples(eventsToAdjust);

dioKinDiffs = corrected_eventSamples./sampling_rate - KinArm_eventTimes;

figure
scatter(KinArm_eventTimes(~missedDigEvents), dioKinDiffs(~missedDigEvents))
hold on
scatter(KinArm_eventTimes(missedDigEvents), dioKinDiffs(missedDigEvents))
legend('Present', 'Missed', 'location', 'Southeast')
title('Kinarm times vs. DIO times (after correction)')

%
disp(['Time drift b/w KinArm and Intan: ', num2str(dioKinDiffs(end)-dioKinDiffs(1)) 's'])
disp(['Total recording time: ' num2str(KinArm_eventTimes(end)) 's'])
disp(['Median lag b/w DIO and UDP: ' num2str(best_lag) ' samples'])

% test_matrix = [KinArm_eventCodes, correctedintcodes,KinArm_eventTimes, corrected_textEventSamples, corrected_eventSamples, corrected_textEventSamples-corrected_eventSamples];

EventTimeSamples = corrected_eventSamples;
Events = KinArm_eventCodes;
