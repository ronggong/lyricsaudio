% Flabread
%
% [t_start t_end label] = Floadlab(file_path)
%
% From Ircam

function [t_start,t_end, label] = Flabread(file_path)

    % Default outputs
    t_start = [];
    t_end = [];
    label = '';

    % Definition of tokens to read from each line of file
    pat = '(?<start>\S+)\s+(?<end>\S+)\s+(?<label>\S+)';

    % Open file
    fid = fopen(file_path,'r');

    % Read data
    while ~feof(fid)
        % read line
        line = fgetl(fid);
        % match tokens
        d = regexp(line,pat,'names');

        t_start(end+1) = str2num(d.start); % was Fstr2num
        t_end(end+1) = str2num(d.end);     % was Fstr2num
        label{end+1} = deblank(d.label);
    end


    % Determine if HTK format (integer values of 100ns units)
    if all(round(t_start)==t_start) & all(round(t_end)==t_end)
        % convert to seconds
        t_start = t_start/10000000;
        t_end = t_end/10000000;
    end

    % Close file
    fclose(fid);

return
