function [F0,har_amp_mat,har_freq_mat] = F0_estimation( ag,hop_int,full_S,tau,cents_stft,Fl,cutoff_curve )
%F0_ESTIMATION estimation of F0 with the possible F0
%   ag : all the possible F0
%   tau : number of stft frame
%   har_amp_mat : matrix of harmonics amplitudes row: harmonics col:time
%   har_freq_mat : matrix of harmonics frequencies
%   F0 : cents unit

N_ag = length(ag);
N_har = 12; % harmonics number
F0 = zeros(1,tau);

ag_mat = zeros(N_ag,tau+1); % matrix used to store the agents frenquencies
err = 20; % error of searching local maximum amplitude

for ii = 1:N_ag
    ag_mat(ii,ag{ii}.time) = ag{ii}.freq*hop_int;
    ag_mat(ii,tau+1) = sum(ag{ii}.reli)/length(ag{ii}.reli);
end

har_amp_mat = zeros(N_har,tau);
har_freq_mat = NaN(N_har,tau);

for ii = 1:tau
    if any(ag_mat(:,ii))
        [r,~,v] = find(ag_mat(:,ii)); % find vertical non zeros frequency
        v_freq = cents2freq(v+Fl-1);
        L = floor(cutoff_curve(3)*(1+err/1200)./v_freq);
        L = min(L);
        N_v = length(v);
        
        amp_mat = zeros(L,N_v); % temp matrix
        freq_mat = zeros(L,N_v); % temp matrix
        for jj = 1:N_v
            for kk = 1:L
                
                [~,low_bound] = min(abs(cents_stft - freq2cents(kk*v_freq(jj))+err));
                [~,up_bound] = min(abs(cents_stft - freq2cents(kk*v_freq(jj))-err));
                
                [amp_mat(kk,jj),freq_mat(kk,jj)] = max(full_S(low_bound:up_bound,ii));
                freq_mat(kk,jj) = freq_mat(kk,jj) + low_bound - 1;
            end
        end
        
        [~,ind_maxpower] = max(sum(amp_mat(1:end,:).^2).*ag_mat(r,tau+1)');
        F0(ii) = v(ind_maxpower);
        
        if  max(L,N_har) == L
            
            har_amp_mat(:,ii) = amp_mat(1:N_har,ind_maxpower);
            har_freq_mat(:,ii) = cents2freq(cents_stft(freq_mat(1:N_har,ind_maxpower)));
        else
            har_amp_mat(1:L,ii) = amp_mat(:,ind_maxpower);
            har_freq_mat(1:L,ii) = cents2freq(cents_stft(freq_mat(:,ind_maxpower)));
        end
    end
end

F0 = F0 + Fl - 1;

end

