function [ line_group ] = prelines( line_group,ii )
%PRELINES preprocess line segments, eliminate overlap line segments in
%                                   groups
%   Detailed explanation goes here
ind_ii = find(line_group(:,5) == ii);
    if length(ind_ii) >=2
        for jj = 1:length(ind_ii)-1
            for kk = jj+1:length(ind_ii)
                temp_line = line_group(ind_ii(jj),:);
                if (line_group(ind_ii(kk),1) <= temp_line(2) && line_group(ind_ii(kk),1)>temp_line(1) && line_group(ind_ii(kk),2)>temp_line(2))
                    
                    line_group(ind_ii(kk),1) =  temp_line(1);
                    line_group(ind_ii(kk),4) = (line_group(ind_ii(kk),4)*line_group(ind_ii(kk),3)+temp_line(4)*temp_line(3))/...
                        (line_group(ind_ii(kk),3)+temp_line(3));
                    line_group(ind_ii(kk),3) =  line_group(ind_ii(kk),2) - temp_line(1);
                    
                    if line_group(ind_ii(kk),2) == line_group(ind_ii(kk),7)          % generative line segment
                        line_group(ind_ii,6) = temp_line(1);
                        line_group(ind_ii,7) = line_group(ind_ii(kk),7);
                    end
                    line_group(ind_ii(jj),:) = nan;
                elseif  (line_group(ind_ii(kk),2) >= temp_line(1) && line_group(ind_ii(kk),1)<temp_line(1) && line_group(ind_ii(kk),2)<temp_line(2))
                    
                    line_group(ind_ii(kk),2) =  temp_line(2);
                    line_group(ind_ii(kk),4) = (line_group(ind_ii(kk),4)*line_group(ind_ii(kk),3)+temp_line(4)*temp_line(3))/...
                        (line_group(ind_ii(kk),3)+temp_line(3));
                    line_group(ind_ii(kk),3) =  temp_line(2) - line_group(ind_ii(kk),1) ;
                    if line_group(ind_ii(kk),1) == line_group(ind_ii(kk),6)          % generative line segment
                        line_group(ind_ii,7) = temp_line(2);
                        line_group(ind_ii,6) = line_group(ind_ii(kk),6);
                    end
                    line_group(ind_ii(jj),:) = nan;
                end
            end
        end
    end

line_group(isnan(sum(line_group,2)),:) = [];
end
