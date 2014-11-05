function [ time,str ] = read_audiosulptmarker( filename )
%READ_MARKER read the label 
%   start_f : event start frame
%   end_f : event end frame
fileID = fopen(filename);
rows = textscan(fileID,'%f %s %s %s %*[^\n]');
time = rows{1}; % in seconds
time = time(:);
%time = cell2mat(time);

str = rows{4};
str = str(:);
fclose(fileID);

end

