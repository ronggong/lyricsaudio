function [ F0_pdf] = EM_KI2( Signal,N,tau,G_W,c_0_1,c_0_2,M,H )
%EM_KI EM algo
%   S: spectrum at time t
%   c_0_1, c_0_2: prior dist amplitude harmonics m = 1,2
%   w_p_1, w_p_2: weight of tone model m = 1,2
% frequencies high and low in cents
% N = length(cents_stft);

hh = waitbar(0,'Estimating melody','Name','Estimating the F0...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        
F0_pdf = zeros(tau,M);
max_it = 10;
loglik_threshold = 1e-3; % converge threshold

c_1_i = zeros(N,M,H);
c_2_i = zeros(N,M,H);
for h = 1:H
    c_1_i(:,:,h) = repmat(c_0_1(h),N,M);
    c_2_i(:,:,h) = repmat(c_0_2(h),N,M);
end

% w : M*N initialise w
w_1 = repmat(1/M/2,1,M);
w_2 = w_1;

    
for t = 1:tau
    %fprintf('t: %i\n',t);
    %tic
    if mod(t,20) == 0
        waitbar(t/tau,hh)
    end
    
    S = Signal(:,t);
    S_NM = repmat(S,1,M);
    S_NMH = repmat(S_NM,[1,1,H]);
    
    % re-initialise c_1
    c_1_N = c_1_i;
    c_2_N = c_2_i;
        
    loglik_old = -realmax;
    converged = 0;
    it = 0;
    
    while ~converged && any(S) % iteration time
        it = it + 1;
        
        %G_W = zeros(M,N,H); % eq 13
        
        ph_1 = c_1_N.*G_W; % eq 11
        ph_2 = c_2_N.*G_W;
%         if any(isnan(ph_1))
%             disp('ph')
%             delete(hh);
%             return;
%         end
        
        %for h = 1:H
        %G_W = zeros(M,N); % eq 13
        %for F = 1:M
        %xx0 = (F+Fl-1)+1200*log2(h);
        
        %G_W(F,:) = exp(-(cents_stft-xx0).^2/2/W^2)/sqrt(2*pi*W^2); % eq 13
        
        % ph_1(F,:,h) = c_1(h,F)*G_W(F,:,h); % eq 11
        % ph_2(F,:,h) = c_2(h,F)*G_W(F,:,h);
        
        %end
        %end
        
        % p N*M
        p_1 = sum(ph_1,3); % eq,10
        p_2 = sum(ph_2,3);
        
        %wp_1 = w_1.*p_1;
        %wp_2 = w_2.*p_2;
        
        wp_1 = bsxfun(@times,p_1,w_1);
        wp_2 = bsxfun(@times,p_2,w_2);
        
%         if ~any(any(wp_1))
%             disp('wp')
%             delete(hh);
%             return;
%         end
%         
        %Norm_factor : N*1
        Norm_factor = sum(wp_1 + wp_2,2); % numerator of eq 43, 44
        Norm_factor(Norm_factor == 0) = eps;
        Nf_NM = repmat(Norm_factor,1,M);
        Nf_NMH = repmat(Nf_NM,[1,1,H]);
%         if ~any(Norm_factor)
%             disp('Norm')
%             delete(hh);
%             return;
%         end
        
        % w_ML : 1*M
        w_ML_1 = sum(S_NM.*(wp_1./Nf_NM)); % eq 43
        w_ML_2 = sum(S_NM.*(wp_2./Nf_NM));
        
        % w : N*M
        w_1 = w_ML_1;
        w_2 = w_ML_2;
        
%         if ~any(any(w_1))
%             disp('w_1')
%             delete(hh);
%             return;
%         end
        
        % c_ML : M*H, eq 44
        c_ML_1 = squeeze(sum(S_NMH.*bsxfun(@times,ph_1,w_1)./Nf_NMH));%./repmat(w_ML_1,H,1);
        c_ML_2 = squeeze(sum(S_NMH.*bsxfun(@times,ph_2,w_2)./Nf_NMH));%./repmat(w_ML_2,H,1);
       
%         if ~any(any(c_ML_1))
%             disp('c_ML')
%             delete(hh);
%             return;
%         end
        %c_1 = (c_ML_1)./repmat(w_ML_1,1,H);
        %c_2 = (c_ML_2)./repmat(w_ML_2,1,H);
        
        c_1 = bsxfun(@rdivide,c_ML_1',w_ML_1);
        c_2 = bsxfun(@rdivide,c_ML_2',w_ML_2);
        
        
        %normalise H*M
        c_1 = bsxfun(@rdivide,c_1,sum(c_1)); % faster than repmat
        c_2 = bsxfun(@rdivide,c_2,sum(c_2));
        %c_1 = c_1./repmat(sum(c_1,2),1,H);
        %c_2 = c_2./repmat(sum(c_2,2),1,H);
        
        % 0/0 = nan, nan is due to 0/(sum = 0)
        c_1(isnan(c_1)) = 0;
        c_2(isnan(c_2)) = 0;
        
        for h = 1:H
            c_1_N(:,:,h) = repmat(c_1(h,:),N,1); % N*M*H
            c_2_N(:,:,h) = repmat(c_2(h,:),N,1);
        end
        % for debug
        %if any(any(isnan(w_1))) || any(any(isnan(w_2))) || any(any(isnan(c_1))) || any(any(isnan(c_2)))
        %    break
        %end
        
        loglik = log(sum(w_1));
        if (abs((loglik/loglik_old)-1) < loglik_threshold) || it > max_it 
            converged = 1;
        end
        
        %fprintf('loglik: %f\n',abs((loglik/loglik_old)-1));
        loglik_old = loglik;
        
    end
    %toc
    %java.lang.Runtime.getRuntime.freeMemory
    %fprintf('t: %i\n',t);
    
    if any(S)
        F0_pdf(t,:) = w_1+w_2; % eq 20
    else
        F0_pdf(t,:) = 0;
    end
end
F0_pdf = F0_pdf';

delete(hh)

end

