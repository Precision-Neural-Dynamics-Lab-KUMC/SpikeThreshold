function [monkey, date_str, arrayString, arrayName, task, file_type] = parse_data_file(fileName)
underscores_i = regexpi(fileName, '_');
period_i = regexpi(fileName, '\.');
if length(underscores_i) >= 2
    monkey = fileName(1:(underscores_i(1)-1));
    date_str = fileName((underscores_i(1)+1):(underscores_i(2)-1));
    arrayString = fileName((underscores_i(2)+1):(underscores_i(3)-1));
    task = fileName((underscores_i(3)+1):(period_i-1));
    file_type = fileName((period_i+1):end);
else
    monkey = [];
    date_str = [];
    arrayString = [];
    task = [];
    file_type = [];
end
non_letters = regexp(arrayString, '[^a-zA-Z]');
arrayName = cellstr(arrayString')';
for k = non_letters
    arrayName{k-1} = [arrayName{k-1} arrayName{k}];
end
arrayName(non_letters) = [];