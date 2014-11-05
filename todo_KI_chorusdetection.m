% this algorithm is based on the article : 
% A Chorus Section Detection Method for Musical Audio Signals and Its Application to a Music Listening Station
% Masataka Goto
% IEEE TRANSACTIONS ON AUDIO, SPEECH, AND LANGUAGE PROCESSING , VOL. 14, NO. 5, SEPTEMBER 2006

% coding by Rong GONG, 5/nov/2014, Paris

%% add the function path 'segmentation'
% load audio file

% addpath(genpath(fullfile(pwd,'..')));
% currentName_AUDIO = fullfile(pwd,'..','04andyliu.wav');
% [aud,fs] = audioread(currentName_AUDIO);

%% fft parameters
fftlen = 0.24*fs; %fftlen = 2048;
hopsize = 0.08*fs; %hopsize = 512; % par defaut, sampling rate 16kHz, 160 samples = 10ms
nFreq = fftlen/2 + 1;

%% Short time fourier transform
[S,freq_stft, time_frame,tau] = STFT_KI( aud,fs,fftlen,hopsize,0 );
S = S(1:nFreq,:);

%% chroma filter and  chromagram
[ chroma_fil ] = chromafil( freq_stft, 48,119 );    % chroma filter
% figure
% imagesc(chroma_fil)
chrogram = zeros(12,tau);
for ii = 1:tau
    chrogram(:,ii) = chroma_fil * S(:,ii);      
end

%% ssm - self-similarity matrix of chromagram
[ ssm_chro ] = chromassm( chrogram );

% imagesc(ssm_chro)
% set(gca,'YDir','normal')
%% remove noise p1786 the last paragraph
[ ssm_chro2,L ] = chromanr( ssm_chro,hopsize,fs );

% imagesc(ssm_chro2)
%% R_all eq.6

R_all_t = sum(ssm_chro2,2)/tau;

%% high pass - removes the global drift caused by cumulative noise
% p1787 paragraph below eq.7

% eliminate boundary abnormal
R_all_t(1:2*L) = R_all_t(2*L);
R_all_t(tau-2*L:tau) = R_all_t(tau-2*L);

B_s = round(16/hopsize*fs);             % 16s window
R_all_t_s = movaverage( B_s,R_all_t );  % moving average
R_all_t_h = R_all_t - R_all_t_s;        % h means high pass

%% find peaks, use True envelop as the threshold

%%% below is the original method in article, eq.7. I didn't use it, because
%%% it's not working
% R_all_t_h_s = zeros(tau,1);
% K_s = round(0.5/hopsize*fs);
% R_all_t_pad = [repmat(R_all_t_h(1),K_s,1);R_all_t_h;repmat(R_all_t_h(end),K_s,1)];
% for ii = 1:tau
%     R_all_t_h_s(ii,1) = sum(R_all_t_pad(ii:ii+2*K_s).*(-K_s:K_s)');
% end

% R_all_diff = diff(R_all_t_h);
% 
% locs = [];
% for ii = 2:tau-1
%     if R_all_diff(ii) < 0 && R_all_diff(ii-1) >0
%         locs(end+1) = ii;
%     end
% end
% 
% % figure
% % plot(1:tau, R_all_t_h,locs,R_all_t_h(locs),'ro')
%%%%%% 

% True envelop parameters
ccOrder = 22050/(2*440);
scaleFactor = 1.3;                              % value can be chosen between 0.5 1.5

R_all_te = [R_all_t_h(1:500,1);R_all_t_h];      % concatenate 500 points in front, because True envelope algorithm can't deal with low frequency range
R_all_t_l = Ftrueenv(abs(R_all_te), ccOrder*scaleFactor);   % l means low pass
R_all_t_l = real(R_all_t_l(501:end));

% figure
% plot(1:tau,R_all_t_h, 1:tau,R_all_t_l)

[peaks,locs] = findpeaks(R_all_t_h,'MINPEAKDISTANCE',2*L,'SORTSTR','descend');      % locs will be flush
peaks_l = R_all_t_l(locs);                                                          % envelope values in peaks location 
locs(peaks<peaks_l | peaks<mean(R_all_t_l)) = [];                                   % eliminate peaks which below threshold and mean value of envelope
peaks(peaks<peaks_l | peaks<mean(R_all_t_l)) = [];
num_peaks = length(locs);

% figure
% plot(1:tau, R_all_t_h,locs,peaks,'ro',1:tau,R_all_t_l)
% legend('R_{all}','peaks','True envelope (threshold)')
% pause
% close gcf


%% search line segments, paragraph of eq.8
th_line = floor(6.4/hopsize*fs);        % threshold of minimum line segment
B_s_p = round(2/hopsize*fs);            % 2s window
line_seg = [];
for ii = 1:num_peaks
    lsf = ssm_chro2(locs(ii),locs(ii):tau)';    % search of line segments will perform on enhanced ssm
    lsf_orig = ssm_chro(locs(ii),locs(ii):tau)';
    lsf_s = movaverage( B_s_p,lsf );        % smoothing
    th = dynthres( lsf_s );                 % dynamic threshold eq.8
    [line_s,line_e,lambda,sigma] = linesegfinder( lsf_s,lsf_orig,th,th_line );    % find line segments
    
    if ~isempty(line_s)
        % line segment matrix, 
        % #1 start time #2 end time #3 lag time #4 mean probability
        % #5 sum prob #6 label 1 - first time detected lines segment (in contrast to the missing lines which is 0)
        line_seg = [line_seg ; [[line_s' line_e']+locs(ii)-1 repmat(locs(ii)-1,length(line_s),1) lambda' sigma', ones(length(line_s),1)]];
    end
    
%     figure
%     plot(1:length(lsf),lsf_s,line_s,lsf_s(line_s),'ro',line_e,lsf_s(line_e),'bo')
%     pause
%     close gcf
end

%%
% figure
% hold on
% for ii = 1:size(line_seg,1)
%     line([line_seg(ii,1) line_seg(ii,2)],[line_seg(ii,3) line_seg(ii,3)])
% end
% hold off

%% III D regroup (integrate) the line segments

T_r = 0.2; T_s = round(3.6/hopsize*fs);
line_group = linesegrouper( ssm_chro, B_s_p,line_seg,T_r,T_s );
l_lg = length(line_group);      % number of groupes

%% 
% the lines with same color are in group
% co = distinguishable_colors(l_lg);
% figure
% hold on
% for ii = 1:l_lg
%     for jj = 1:size(line_group{ii,1})
%         line([line_group{ii,1}(jj,1) line_group{ii,1}(jj,2)],[line_group{ii,1}(jj,3) line_group{ii,1}(jj,3)],'color',co(ii,:))
%     end
% end
% hold off

%% redetect missing segment eq.9, same method to finding line segments
R_T = zeros(tau,l_lg);
R_T_h = R_T;
for ii = 1:l_lg
    G_s = line_group{ii,2}(1);      % group start 
    G_e = line_group{ii,2}(2);      % group end
    
    % for the raison of find peaks, the denominator is modified to tau
    % instead of G_e-G_s+1 in eq.9
    R_T(:,ii) = sum(ssm_chro2(:,G_s:G_e),2)/tau;
    R_T(1:2*L,ii) = R_T(2*L,ii);                    % eliminate boundary abnormal
    [ R_T_s ] = movaverage( B_s,R_T(:,ii) );
    R_T_h(:,ii) = R_T(:,ii) - R_T_s;
end

%% find peaks use low pass envelope as the threshold
locs_R_T = cell(l_lg,1);
for ii = 1: l_lg
[peaks,locs] = findpeaks(R_T_h(:,ii),'MINPEAKDISTANCE',L,'SORTSTR','descend');

B_s_th = round(4/hopsize*fs);                   % 4s window
R_T_l = movaverage( B_s_th,R_T_h(:,ii) );       % threshold
peaks_l = R_T_l(locs)+0.005;
locs(peaks<peaks_l) = [];
peaks(peaks<peaks_l) = [];
locs_R_T{ii,1} = [peaks locs];                  % peaks vectors
            
%     figure
%     plot(1:tau, R_T_h(:,ii),locs, peaks,'ro')
%     pause
%     close gcf

end
%% integrate missing line segments into group p1778
U_r = 1.4;
for jj = 1:l_lg
    G_s = line_group{jj,2}(1);
    G_e = line_group{jj,2}(2);
    if ~isempty(locs_R_T{jj,1})
        line_seg = [];
        temp = locs_R_T{jj,1};
        
        % removes inappropriate peaks Rule 2), standard deviation threshold
        temp_smm = ssm_chro(line_group{jj,1}(:,3),G_s:G_e);
        for ii = 1:size(temp_smm,1)
            temp_smm(ii,:) = movaverage( B_s_p,temp_smm(ii,:)' )';
        end
        temp(std(temp_smm,0,2)>U_r*line_group{jj,3},:) = [];
        %%%
        
        for ii = 1:size(temp,1)
            lsf = ssm_chro2(temp(ii,2),G_s:G_e)';
            lsf_orig = ssm_chro(temp(ii,2),G_s:G_e)';
            [ lsf_s ] = movaverage( B_s_p,lsf );
            [ th ] = dynthres( lsf_s );
            if ~isempty(th)
                [ line_s,line_e,lambda,sigma ] = linesegfinder( lsf_s,lsf_orig,th,th_line );
            end
%                 figure
%                 plot(1:length(lsf),lsf_s,line_s,lsf_s(line_s),'ro',line_e,lsf_s(line_e),'bo')
%                 pause
%                 close gcf
            if ~isempty(line_s)
                line_seg = [line_seg ; [[line_s' line_e']+G_s-1 repmat(temp(ii,2)-1,length(line_s),1) lambda' sigma' 2*ones(length(line_s),1)]];
            end
        end
        
        for kk = 1: size(line_seg,1)
            temp_line = line_seg(kk,:);
            if ~any(sum(abs(bsxfun(@minus,line_group{jj,1},temp_line))) == 0)   % jugde if line_seg already in line group
                line_group{jj,1} = [line_group{jj,1};temp_line];                % add missing line segments to group
            end
        end
    end
end


%% removes inappropriate peaks (missing line segments) Rule 3) 
%  Remove a peak that is too close to other peaks and causes sections to overlap.

for ii = 1: l_lg
    G_s = line_group{ii,2}(1);
    G_e = line_group{ii,2}(2);
    
    [~,ind_order ] = sort(line_group{ii,1}(:,5),'descend');     % sort line segments in each group regarding their sum probabilities 
    line_group{ii,1} = line_group{ii,1}(ind_order,:);
   
    jj = 1;
    while 1
        rel_dis = abs(line_group{ii,1}(:,3) - line_group{ii,1}(jj,3));  % distance between two peaks
        ind_rem = find(rel_dis<G_e-G_s+1 & line_group{ii,1}(:,6) ~= 1); % not remove the line segments of first detection
        ind_rem(ind_rem == jj) = [];
        line_group{ii,1}(ind_rem,:) = [];
        
        if jj >= size(line_group{ii,1},1)
            break
        end
        jj = jj + 1;
    end
end


%%
% figure
% hold on
% for ii = 1:l_lg
%     for jj = 1:size(line_group{ii,1})
%         line([line_group{ii,1}(jj,1) line_group{ii,1}(jj,2)],[line_group{ii,1}(jj,3) line_group{ii,1}(jj,3)],'color',co(ii,:))
%     end
% end
% hold off

%% unfold line groups to a matrix
kk = 0;
line_group_unfold_mat = [];
for ii = 1:l_lg
    G_s = line_group{ii,2}(1);
    G_e = line_group{ii,2}(2);
    
    for jj = 1:size(line_group{ii,1})
        kk = kk + 1;
        temp_s = line_group{ii,1}(jj,1);    % start time
        temp_e = line_group{ii,1}(jj,2);    % end time
        temp_l = line_group{ii,1}(jj,3);    % length
        
        % assumption 1 p1789
        if temp_e-temp_s > 7.7/hopsize*fs && temp_e-temp_s < 50/hopsize*fs
            temp_lam = line_group{ii,1}(jj,4);
        else
            temp_lam = 0;
        end
       
        % #1 start #2 end #3 length #4 lambda (mean probability) #5 group number
        % #6 group start bin #7 group end bin #8 label of concatenation, 0 - not yet concatenated
        
        line_group_unfold_mat(end+1,:) = [temp_s - temp_l, temp_e - temp_l, temp_e - temp_s,temp_lam, ii, G_s, G_e,0];
    end
    
    % generative line segment
    if G_e - G_s > 7.7/hopsize*fs && G_e - G_s < 50/hopsize*fs
        max_lam = max(line_group{ii,1}(:,4));
    else
        max_lam = 0;
    end
    line_group_unfold_mat(end+1,:) = [G_s, G_e, G_e-G_s,max_lam, ii, G_s,G_e,0];
    
end

% eliminate the line segments which mean probabilities are < th_prob
th_prob = 0.5;
line_group_unfold_mat(line_group_unfold_mat(:,4)<th_prob,:) = [];


%% eliminate overlap line segments in groups

for ii = 1:l_lg
    line_group_unfold_mat = prelines( line_group_unfold_mat,ii );
end
%% regroup p1788 below rule 3)

line_group_unfold_mat = linesregrouper( line_group_unfold_mat, l_lg, hopsize, fs );

%% post post processing
% remove short line segments < 1/3 original group length
for ii = 1:l_lg
    ind_ii = find(line_group_unfold_mat(:,5) == ii);
    
    if length(ind_ii) >= 2
        G_s = line_group{ii,2}(1);
        G_e = line_group{ii,2}(2);
        
        ind_rem = line_group_unfold_mat(ind_ii,3) < (G_e-G_s)/3;
        
        line_group_unfold_mat(ind_ii(ind_rem),:) = [];
    end
    
end

% remove inner lines in groups

for ii = 1:l_lg
    ind_ii = find(line_group_unfold_mat(:,5) == ii);
    
    if length(ind_ii) > 2
        ind_rem = [];
        
        for jj = 1:length(ind_ii)
            ind_temp = ind_ii(line_group_unfold_mat(ind_ii,1) >= line_group_unfold_mat(ind_ii(jj),1) & line_group_unfold_mat(ind_ii,2) <= line_group_unfold_mat(ind_ii(jj),2));
            ind_temp(ind_temp == ind_ii(jj)) = [];
            ind_rem = [ind_rem; ind_temp];
        end
        line_group_unfold_mat(ind_rem,:) = [];
    end
    
end

%% assumption 2 p1789

T_long = 50/hopsize*fs;
[max_seg_l, max_seg_i] = max(line_group_unfold_mat(:,3));
upper_num = sum(line_group_unfold_mat(:,3) > T_long);
if upper_num
    upper_ind = find(line_group_unfold_mat(:,3) > T_long);                              % find non zeros indices
    temp_unfold_mat = line_group_unfold_mat(line_group_unfold_mat(:,3) <= T_long,:);    % segments length < T_long
    temp_unfold_mat_ind = find(line_group_unfold_mat(:,3) <= T_long);
    for ii = 1:length(upper_ind)
        temp_ind = upper_ind(ii);
        temp_end_base = line_group_unfold_mat(temp_ind,2);                              % end time to compare
        for jj = 1:size(temp_unfold_mat,1)
            if abs(temp_unfold_mat(jj,2) - temp_end_base) < T_s
                line_group_unfold_mat(temp_unfold_mat_ind(jj),4) = temp_unfold_mat(jj,4)*2;
            end
        end
    end
end

%% assumption 3 p1789
l_lgu = size(line_group_unfold_mat,1);
for ii = 1:l_lgu
    line_seg_base = line_group_unfold_mat(ii,:);        % line for compare
    
    temp_l = line_group_unfold_mat(:,3);
    temp_s = line_group_unfold_mat(:,1);
    temp_e = line_group_unfold_mat(:,2);
    
    % find a segment is halflength, and between the section
    temp_ind = find(temp_l < line_seg_base(3)/2 & temp_s > line_seg_base(1) & temp_e < line_seg_base(2));
    
    if length(temp_ind) >= 2
        lu = unique(line_group_unfold_mat(temp_ind,5));
        if length(lu) == 1, counts = length(line_group_unfold_mat(temp_ind,5));
        else counts = hist(line_group_unfold_mat(temp_ind,5),lu); end
        
        if any(counts == 2)
            ind_rep = lu(counts == 2);                  % group number which repeat 2 times
            ind_jj = [];
            for jj = 1:length(ind_rep)
                ind_jj = [ind_jj;find(line_group_unfold_mat(temp_ind,5) == ind_rep(jj))];
                
            end
            
            line_group_unfold_mat(ii,4) = line_group_unfold_mat(ii,4) + mean(line_group_unfold_mat(temp_ind(ind_jj),4))/2;
        end
    end
end
%% select chorus sections F
v = zeros(l_lg,1);
D_l = 1.4/hopsize*fs;
for ii = 1:l_lg
    ind_ii = find(line_group_unfold_mat(:,5)== ii);
    if ~isempty(ind_ii)
        v(ii) = sum(line_group_unfold_mat(ind_ii,4).*line_group_unfold_mat(ind_ii,3))*log((mean(line_group_unfold_mat(ind_ii,3))+1)/D_l);
%         v(ii) = sum(line_group_unfold_mat(ind_ii,4))*log((mean(line_group_unfold_mat(ind_ii,3))+1)/D_l);
    else
        v(ii) = 0;
    end
end
[~,m_can] = sort(v,'descend');      % chorus candidats
m = m_can(1);                       % chorus group

% chorus usually is the last section
% if max(line_group_unfold_mat(line_group_unfold_mat(:,5) == m_can(1),2)) > ...
%     max(line_group_unfold_mat(line_group_unfold_mat(:,5) == m_can(2),2))
%     m = m_can(1);
% else
%     m = m_can(2);
% end
%%
figure
hold on
for jj = 1:l_lg
    ind_m = find(line_group_unfold_mat(:,5)== jj);

    if ~isempty(ind_m)
        for ii = 1:length(ind_m)
            temp_x1 = time_frame(line_group_unfold_mat(ind_m(ii),1));
            temp_x2 = time_frame(line_group_unfold_mat(ind_m(ii),2));
            line([temp_x1 temp_x2],[jj jj]);
        end
    end

end
ylim([0 l_lg+1])
% vline(time_fri,repmat({'r:'},length(time_fri),1),str_fri) % truth
title(currentName_AUDIO)
hold off

chorus_section = line_group_unfold_mat(line_group_unfold_mat(:,5) == m,1:2);
chorus_section = time_frame(chorus_section);

