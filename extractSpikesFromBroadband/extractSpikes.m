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


function version = extractSpikes(dataPaths, envInfo, dataBlocks, filtInfo)

if ischar(dataPaths)
    version = '2.0';
    return
end

if nargin < 4 || isempty(filtInfo)
    filtInfo.filt_order = 4;  %4th order filter
    filtInfo.band_limits = [250, 5000]; % bandpass between 250-7500 Hz
    filtInfo.time_pre       = 175;    % Amount of time before trigger for snippet (microseconds)
    filtInfo.time_post      = 625;    % Amount of time after trigger for snippet (microseconds)
    filtInfo.time_peak_excl = 625;   %Minimum time from previous threshold crossing that the next spike can occur
    filtInfo.time_req_baseline = 175;  %Minimum time signal must be below threshold crossing before next spike can occur
    filtInfo.peak_window    = 150;   %Time after trigger where waveform peak can occur  (microseconds)
    filtInfo.align_spikes   = false;
    filtInfo.throwout_crosstalk = false;
    filtInfo.throwout_large_artifact = false;
end

if dataPaths.intan
    if dataPath.RHD
        extractSpikesRHD(dataPaths, envInfo, dataBlocks, filtInfo)
    else
        extractSpikesRHS(dataPaths, envInfo, dataBlocks, filtInfo)
    end
end
if Datafiles.spikeGadgets
    extractSpikesREC(dataPaths, envInfo, dataBlocks, filtInfo)
end
if ~(dataPaths.intan  || Datafiles.spikeGadgets)
    extractSpikesNsX(dataPaths, envInfo, dataBlocks, filtInfo)  %Note, for backwards compatability but this has not been tested thoroughly
end
