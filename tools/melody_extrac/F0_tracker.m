function ag = F0_tracker(F0_pdf,M)
% ch 3.3
% F0_pdf: row, Freq; colomn: time

min_peak_dis = floor(M/10);
[p,locs] = findpeaks(F0_pdf(:,1),'MINPEAKDISTANCE',min_peak_dis,'SORTSTR','descend');

if length(p) > 3
    p = p(1:3);
    locs = locs(1:3);
end

% initialise agents
th = p(1)/2; % threshold (most dominant peak/2)
N_agent = length(find(p > th)); % amont of agent

for ii = 1:N_agent
    ag{ii}.time = 1;
    ag{ii}.freq = locs(ii);
    ag{ii}.reli = p(ii)/p(1);
    ag{ii}.pena = 0;
end

ag_alive = (1:N_agent); % live agent
last_ag = N_agent; % the number of last agent
ag_penalty = []; % penalty agent, vide intial

for t = 2:size(F0_pdf,2)
    
    [p,locs] = findpeaks(F0_pdf(:,t),'MINPEAKDISTANCE',min_peak_dis,'SORTSTR','descend');
    
    % search peak directly in PDF, step 3)
    if ~isempty(ag_penalty)
       for ii = 1:length(ag_penalty)
           d = abs(locs - ag{ag_penalty(ii)}.freq(end));
           [min_d,J] = min(d); % find the min distance between agent and peaks
           % search the peak in the interval of 200 points
           if min_d < min_peak_dis/2 % if found
               ag{ag_penalty(ii)}.time = [ag{ag_penalty(ii)}.time t];
               ag{ag_penalty(ii)}.freq = [ag{ag_penalty(ii)}.freq locs(J)];
               ag{ag_penalty(ii)}.reli = [ag{ag_penalty(ii)}.reli p(J)/p(1)];
               ag{ag_penalty(ii)}.pena = 0;
               ag_penalty(ii) = 0;
           else % if not found peak, penalize
               ag{ag_penalty(ii)}.pena = ag{ag_penalty(ii)}.pena + 1;
               if ag{ag_penalty(ii)}.pena > 3
                   ag_alive(ismember(ag_alive,ag_penalty(ii)))=[];
               end
           end
       end
       ag_penalty(ag_penalty==0) = []; % if found, delete it from the penalty vector
    end
    
    % limit the peak number to 3 
    if length(p) > 3
        p = p(1:3);
        locs = locs(1:3);
    end
   
    th = p(1)/2;
    N_peak = length(find(p > th));
    
    % 
    ag_non_penalty = ag_alive;
    ag_non_penalty(ismember(ag_alive,ag_penalty)) = []; % agent non penalty = agent alive - penalty
    
    ag_penalty = ag_alive; % re-intialise the agents penalty, if agent attribute after, exclude from this list
    
    peak_alloc = []; % for count which peak is allocated
    I_previous = []; % for count which agent is attributed
    
    % calculate the frequency closeness, step 2)
    d = zeros(length(ag_non_penalty),N_peak);
    for ii = 1:length(ag_non_penalty)
        for jj = 1:N_peak
            d(ii,jj) = abs(ag{ag_non_penalty(ii)}.freq(end)-locs(jj));
        end
    end
    
    % attribute peaks to agents step 2)
    % rule: attribute always the most close first,
    % if an agent is alreay attributed, it can't be attributed again
    % 
    for ii = 1:min(N_peak,length(ag_non_penalty))
        [~,min_ind] = min(d(:));
        [I,J] = ind2sub(size(d),min_ind);
 
        if ~any(ismember(I_previous,I))
            ag{ag_non_penalty(I)}.time = [ag{ag_non_penalty(I)}.time t];
            ag{ag_non_penalty(I)}.freq = [ag{ag_non_penalty(I)}.freq locs(J)];
            ag{ag_non_penalty(I)}.reli = [ag{ag_non_penalty(I)}.reli p(J)/p(1)];
            d(:,J) = Inf;
            
            % attributed agent, exclude from penalty list
            ag_penalty(ismember(ag_penalty,ag_non_penalty(I)))=[];
            peak_alloc = [peak_alloc J]; % peak allocated, add to allocated list
        else
            d(:,J) = Inf;
        end
        I_previous = [I_previous I]; % agent attributed, add to attributed list
    end
    
    if ~isempty(ag_penalty)
        for ii = 1:length(ag_penalty)
            ag{ag_penalty(ii)}.pena = ag{ag_penalty(ii)}.pena + 1;
            %if ag{ag_penalty(ii)}.pena > 3
            %    ag_alive(ismember(ag_alive,ag_penalty(ii)))=[];
            %end
        end
    end
    
    % generate new agent step 2)
    if ~any(ismember(peak_alloc,1))
        ag{last_ag+1}.time = t;
        ag{last_ag+1}.freq = locs(1);
        ag{last_ag+1}.reli = 1;
        ag{last_ag+1}.pena = 0;
        ag_alive = [ag_alive last_ag+1];
        last_ag = last_ag + 1;
    end
    
end

end