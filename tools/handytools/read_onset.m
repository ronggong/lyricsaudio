function [ frames ] = read_onset( filename,time_frame )
%READ_MARKER read the onset marker 
%   
%   
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


end

