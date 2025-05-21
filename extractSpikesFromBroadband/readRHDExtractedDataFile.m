function data_strut = readRHDExtractedDataFile(filename,t_start,t_end,varargin)

% read "One File Per Channel" RHD file
[amplifier_channels, frequency_parameters] = read_Intan_RHD2000_file(dir(filename).folder);
% num_channels = length(amplifier_channels); % amplifier channel info from header file
fileinfo = dir(filename);
num_samples = fileinfo.bytes/(2); % int16 = 2 bytes
fid = fopen(filename, 'r');
%switch variable inputs
switch nargin 
    case 1
        tmp_data = fread(fid, num_samples, 'int16');
    case 3
        fseek(fid,(t_start-1)*2,'cof');
        tmp_data = fread(fid, [1,t_end - t_start + 1],'int16');
    otherwise
        error('unexpected inputs, either only ')
end
data_strut.fields.data = tmp_data';
data_strut.clockrate = frequency_parameters.amplifier_sample_rate;

fclose(fid);
