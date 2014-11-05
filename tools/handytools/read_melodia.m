function [ time,freq ] = read_melodia( filename )
%READ_MELODIA for reading melodia markers
%   time and freq

fileID = fopen(filename);
rows = textscan(fileID,'%f %f %*[^\n]');
time = rows{1}; % in seconds
time = time(:);
freq = rows{2};
freq = freq(:);

fclose(fileID);

end

