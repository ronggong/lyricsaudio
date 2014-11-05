function [ line_group ] = linesegrouper( ssm_chro,B_s_p ,line_seg,T_r,T_s )
%LINESEGROUPER regroupe the line segment
%   goto 2006 III section D
% input:
%   ssm_chro : self similarity matrix chroma
%   B_s_p :    window for moving avearge
%   line_seg : line segments
%   T_r, T_s : T_ratio, T_size, the threshold for combine line segments
%              into group
% output:
%   line_group{:,1} : line segments 
%   line_group{:,2} : group start time, group end time
%   line_group{:,3} : probabiltiy standard deviation
line_group = cell(size(line_seg,1),3);
ii = 1;
while ~isempty(line_seg)
    line_base = line_seg(1,:);              % take first line as base line to compare with
    l_lb = line_base(2) - line_base(1);     % length of base line
    line_group{ii,1}(end+1,:) = line_base;  % add base line to the group
    line_seg(1,:) = [];
    
    if isempty(line_seg)
        break
    else
        line_ind = [];
        for jj = 1:size(line_seg,1)
            l_lc = line_seg(jj,2) - line_seg(jj,1);     % length of line to compare
            th_T = min(T_r*min(l_lb,l_lc), T_s);
            if abs(line_base(1) - line_seg(jj,1)) < th_T && abs(line_base(2) - line_seg(jj,2) < th_T)
                line_group{ii,1}(end+1,:) = line_seg(jj,:);
                line_ind(end+1) = jj;
            end
        end
        line_seg(line_ind,:) = [];                      % delete all the grouped lines
    end
    ii = ii + 1;
end

line_group = line_group(~cellfun('isempty',line_group));

% add start, end point and standard deviation to each group
for ii = 1:size(line_group,1)
    
    line_group{ii,2} = [min(line_group{ii,1}(:,1)) max(line_group{ii,1}(:,2))];     % start, end
    
    temp = ssm_chro(line_group{ii,1}(:,3),line_group{ii,2}(1):line_group{ii,2}(2));
    for jj = 1:size(temp,1)
        temp(jj,:) = movaverage( B_s_p,temp(jj,:)' )';
    end
    std_base = max(std(temp,0,2));
    
    line_group{ii,3} = std_base;        %standard deviation
    
end


end

