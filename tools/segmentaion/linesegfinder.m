function [ line_s,line_e,lambda,sigma ] = linesegfinder( lsf_s,lsf,th,th_line )
%LINESEGFINDER find the line segments
% input:
%   lsf_s : smoothed input signal
%   lsf:    not smoothed input signal
%   th :    threshold above it is line segment, below no
%   th_line : minimun length of line seg
% output
%   line_s, line_e : start and end point of a segment
%   lambda : mean similarity probabiltiy
%   sigma : sum similarity probability


    line_s = []; % line start
    line_e = [];
    if ~isempty(th)
        
        for jj = 2:length(lsf_s)
            if lsf_s(jj-1) < th && lsf_s(jj) >=th
                line_s(end+1) = jj-1;
            end
            if lsf_s(jj-1) >= th && lsf_s(jj) <th && ~isempty(line_s)
                line_e(end+1) = jj;
            end
        end
        if length(line_s) > length(line_e)
            line_s(end) = [];
        end
    end
   
    ii = 2;
    while 1
        if ii >= length(line_s)
            break
        end
        if line_s(ii) -line_e(ii-1) == 1 && line_e(ii) - line_s(ii) >= 10
            line_e(ii-1) = line_e(ii);
            line_s(ii) = [];
            line_e(ii) = [];
        else
            ii = ii + 1;
        end
    end
    
    if ~isempty(line_s)
        % select the line segments > th_line
        line_diff = line_e - line_s;
        line_s = line_s(line_diff > th_line);
        line_e = line_e(line_diff > th_line);
    else
        line_s = [];
        line_e = [];
    end
    
    % mean probabilty
    if ~isempty(line_s)
    
        lambda = zeros(1,length(line_s));
        sigma = zeros(1,length(line_s));
        for ii = 1:length(line_s)
            lambda(ii) = mean(lsf(line_s(ii):line_e(ii)));
            sigma(ii) = sum(lsf(line_s(ii):line_e(ii)));
        end
    else
        lambda = [];
        sigma = [];
    end
     

end

