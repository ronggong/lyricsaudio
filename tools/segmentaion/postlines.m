function [ line_group,temp_group ] = postlines( line_group,temp_group,ii,th )
%postlines post processing, concatenate overlapping lines segments or which end to start distance < th
             
% line_group, temp_group
% ii:group number

 ind_ii = find(line_group(:,5) == ii);
    if length(ind_ii) >=2
        for jj = 1:length(ind_ii)-1
            for kk = jj+1:length(ind_ii)
                temp_line = line_group(ind_ii(jj),:);
                if (line_group(ind_ii(kk),1) - temp_line(2)  < th && temp_line(2) <= line_group(ind_ii(kk),1))||...
                    (line_group(ind_ii(kk),1) < temp_line(2) && line_group(ind_ii(kk),1)>temp_line(1) && line_group(ind_ii(kk),2)>temp_line(2))
                %debug
%                     temp_line(1)
%                     temp_line(6)
%                     line_group(ind_ii(kk),2)
%                     line_group(ind_ii(kk),7)
%                     disp('1')
%                     pause
                    line_group(ind_ii(kk),1) =  temp_line(1);

                    line_group(ind_ii(kk),4) = (line_group(ind_ii(kk),4)*line_group(ind_ii(kk),3)+temp_line(4)*temp_line(3))/...
                        (line_group(ind_ii(kk),3)+temp_line(3));
                    line_group(ind_ii(kk),3) =  line_group(ind_ii(kk),2) - temp_line(1);
                    
                    
                    if line_group(ind_ii(kk),2) == line_group(ind_ii(kk),7)          % generative line segment

                        line_group(ind_ii,6) = temp_line(1);
                        line_group(ind_ii,7) = line_group(ind_ii(kk),7);
  
                    end
                    line_group(ind_ii(jj),:) = nan;
                    temp_group(ind_ii(jj),:) = nan;
                elseif  temp_line(1) - line_group(ind_ii(kk),2) < th  && temp_line(1) >= line_group(ind_ii(kk),2)||...
                    (line_group(ind_ii(kk),2) > temp_line(1) && line_group(ind_ii(kk),1)<temp_line(1) && line_group(ind_ii(kk),2)<temp_line(2))
% debug                    
%                     line_group(ind_ii(kk),1)
%                     line_group(ind_ii(kk),6)
%                     temp_line(2)
%                     temp_line(7)
%                     disp('2')
%                     pause
                    line_group(ind_ii(kk),2) =  temp_line(2);

                    
                    line_group(ind_ii(kk),4) = (line_group(ind_ii(kk),4)*line_group(ind_ii(kk),3)+temp_line(4)*temp_line(3))/...
                        (line_group(ind_ii(kk),3)+temp_line(3));
                    line_group(ind_ii(kk),3) =  temp_line(2) - line_group(ind_ii(kk),1) ;
                    
                    if line_group(ind_ii(kk),1) == line_group(ind_ii(kk),6) % generative line segment
                        line_group(ind_ii,7) = temp_line(2);
                        line_group(ind_ii,6) = line_group(ind_ii(kk),6);
                    end
                    
                    line_group(ind_ii(jj),:) = nan;
                    temp_group(ind_ii(jj),:) = nan;
                end
            end
        end
    end
    
end


