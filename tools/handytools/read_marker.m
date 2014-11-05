function [ start_f,end_f ] = read_marker( filename,time_frame )
%READ_MARKER read the label 
%   start_f : event start frame
%   end_f : event end frame
fileID = fopen(filename);
rows = textscan(fileID,'%f %*[^\n]');
time = rows{1}; % in seconds
time = (time(:))';
fclose(fileID);

N = length(time);
frames = zeros(1,N);

for ii = 1:N
    [~,frames(ii)] = min(abs(time_frame - time(ii)));
end

start_f = frames(1:2:end);
end_f = frames(2:2:end);


end

