function [ line_group ] = linesregrouper( line_group, l_lg, hopsize, fs )
%LINESREGROUPER regroup p1788 below rule 3)
%   line_group : line segments matrix
%   l_lg : number of groups

T_r_p = 0.1; T_s_p = 3/hopsize*fs; th_T3 = 3.6/hopsize*fs; th_T4 = 23.1/hopsize*fs;
for ii = 1:l_lg
 group_base = line_group(line_group(:,5) == ii,:);  % base group to be compared with
    
    if ~isempty(group_base)
        
        temp_group_mat = line_group;
        
        jj = 0;
        while 1
            if jj >= size(group_base,1)
                break
            end
            jj = jj+1;
            
            % regroup rules : 1) group section is overlapping with the base
            %                 line segment
            %                 2) one group section include anthoer
            %                 3) two group section are close but not
            %                 overlaps. besides, lines segments length are
            %                 < th_T4
            th_T = min(T_r_p*group_base(jj,3), T_s_p);
            
            diff_l = abs(temp_group_mat(:,6) - group_base(jj,1)); diff_l(diff_l==0) = inf;
            diff_r = abs(temp_group_mat(:,7) - group_base(jj,2)); diff_r(diff_r==0) = inf;
            
            diff_l2 = abs(temp_group_mat(:,6) - group_base(jj,6)); diff_l2(diff_l2==0) = inf;
            diff_r2 = abs(temp_group_mat(:,7) - group_base(jj,7)); diff_r2(diff_r2==0) = inf;
            
            diff_r3 = temp_group_mat(:,6) - group_base(jj,7);
            
            ind_jj = find( (diff_l < th_T & diff_r < th_T)...
                | (temp_group_mat(:,6) > group_base(jj,1) &...
                 temp_group_mat(:,7) < group_base(jj,2))...
                 | (diff_l2 < th_T | diff_r2 < th_T));
             
            if group_base(jj,8) == 0                % if jj has been never concatenated
                ind_jj3 = find(diff_r3<th_T3 & diff_r3 > 0 & temp_group_mat(:,3) < th_T4 ...
                    & group_base(jj,3) < th_T4);
            else
                ind_jj3=[];
            end
             
% debug           
%               ii
%               jj
%            diff_l
%               diff_r
%              (line_group_unfold_mat(:,6) > group_base(jj,1) &...
%                  line_group_unfold_mat(:,7) < group_base(jj,2))
%               ind_jj3
%                diff_l2
%                  diff_r2
              
%             
%             temp_group_mat(ind_jj,6)
%             group_base(jj,1)
%                  temp_group_mat(ind_jj,7)
%                  group_base(jj,2)
%              line_group_unfold_mat( ind_jj,5)
%              
%              temp_group_mat(temp_group_mat(:,5) == 29 ,:)
%                pause
%              group_base(jj,7)
%              temp_group_mat(:,6)
             
             

             line_group(ind_jj3,5) = ii;
             line_group(ind_jj,5) = ii;                   % modify the group number
             line_group(ind_jj,6) = group_base(jj,6);     % modify group start
             line_group(ind_jj,7) = group_base(jj,7);     % modify group end
             
             temp_group_mat(ind_jj,:) = nan;
             temp_group_mat(ind_jj3,:) = nan;
           
             if ~isempty(ind_jj3)                                    % if group ii has been concatenated, set concatenation label to 1
                 line_group(line_group(:,5)==ii,8) = 1;
             elseif any(line_group(line_group(:,5)==ii,8))    % if group ii has any concatenation
                 line_group(line_group(:,5)==ii,8) = 1;
             end
             
             % check if there is overlap line
%              group_base = line_group_unfold_mat(line_group_unfold_mat(:,5) == ii,:);
%              temp_line = zeros(max(group_base(:,2)),1);
%              for kk = 1:size(group_base,1)
%                  temp_line(group_base(kk,1):group_base(kk,2)) = 1;
%              end
             
%              sum(temp_line)
%              sum(group_base(:,2)-group_base(:,1)+1)
%              ii
%              jj
%              pause
            
            % concatenate close/overlapping lines in groups
            if ~isempty(ind_jj) || ~isempty(ind_jj3)
                [ line_group,temp_group_mat ] = postlines( line_group,temp_group_mat,ii,th_T3 );
            end
            
            group_base = line_group(line_group(:,5) == ii,:);
            group_base(find(1-any(group_base,2)),:) = nan;
            
            %plot debug
     
%             figure
%             hold on
%             for n = 1:l_lg
%                 ind_m = find(line_group_unfold_mat(:,5)== n);
%                 
%                 if ~isempty(ind_m)
%                     for kk = 1:length(ind_m)
% %                         temp_x1 = time_frame(line_group_unfold_mat(ind_m(kk),1));
% %                         temp_x2 = time_frame(line_group_unfold_mat(ind_m(kk),2));
%                         temp_x1 = line_group_unfold_mat(ind_m(kk),1);
%                          temp_x2 = line_group_unfold_mat(ind_m(kk),2);
%                         line([temp_x1 temp_x2],[n n]);
%                     end
%                 end
%                 
%             end
%             ylim([0 l_lg+1])
% %             vline(time_fri,repmat({'r:'},length(time_fri),1),str_fri) % truth
%             hold off
%             pause
%             if mod(jj,2) == 0
%             close gcf
%             end
        end
    end
end
line_group(isnan(sum(line_group,2)),:) = [];
end

