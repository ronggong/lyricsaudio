function audiosulptlab_convert( filename , writefile, type)
% change the format of audiosulpt lab file
% Default outputs
%t_start = [];

% convert type 1
% example : 3.4444 null null label_1
%           4.3333 null null label_2
%           5.6666 null null label_3

%                      |

%           0       3.4444  label_1
%           3.4444  4.3333  label_2
%           4.3333  5.6666  label_3

% convert type 2
% example : 3.4444 null null label_1
%           4.3333 null null 
%           5.6666 null null label_2
%           4.5555 null null 

%                      |

%           3.4444  4.3333  label_1
%           5.6666  4.5555  label_2


%label = '';

% Open file
fid = fopen(filename,'r');
fid2 = fopen(writefile,'a+');

if type == 1
    t_end = [];
    % Definition of tokens to read from each line of file
    pat = '(?<end>\S+)\s+(?<null1>\S+)\s+(?<null2>\S+)\s+(?<label>\S+)';
    
    ii = 1;
    % Read data
    while ~feof(fid)
        % read line
        line = fgetl(fid);
        % match tokens
        d = regexp(line,pat,'names');
        
        if ii == 1
            fprintf(fid2, '%f %f %s\n', 0, str2num(d.end), deblank(d.label));
            %t_start(end+1) = 0;
            t_end(end+1) = str2num(d.end);
            %label{end+1} = deblank(d.label);
        end
        
        fprintf(fid2, '%f %f %s\n', t_end(end), str2num(d.end), deblank(d.label));
        %t_start(end+1) = t_end(end); % was Fstr2num
        t_end(end+1) = str2num(d.end);     % was Fstr2num
        %label{end+1} = deblank(d.label);
        
        ii = ii + 1;
    end

elseif type == 2
    % Definition of tokens to read from each line of file
    pat1 = '(?<marker>\S+)\s+(?<null1>\S+)\s+(?<null2>\S+)\s+(?<label>\S+)';
    pat2 = '(?<marker>\S+)\s+(?<null1>\S+)\s+(?<null2>\S+)';
    ii = 1;
    % Read data
    while ~feof(fid)
        % read line
        line = fgetl(fid);
        % match tokens
        
        
        if mod(ii,2) == 0
            d2 = regexp(line,pat2,'names');
            
            fprintf(fid2, '%s %s %s\n', d1.marker,d2.marker, 'fric');
        else
            d1 = regexp(line,pat1,'names');
            if ii == 1
                fprintf(fid2, '%f %s %s\n', 0,d1.marker, 'nonfric')
            else
                fprintf(fid2, '%s %s %s\n', d2.marker,d1.marker, 'nonfric') 
            end
           
        end
        
        ii = ii + 1;
    end
    
end
fclose(fid);
fclose(fid2);

end

