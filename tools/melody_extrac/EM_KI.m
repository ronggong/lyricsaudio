function [ F0_pdf] = EM_KI( Signal,N,tau,G_W,c_0_1,c_0_2,M,H )
%EM_KI EM algo
%   S: spectrum at time t
%   c_0_1, c_0_2: prior dist amplitude harmonics m = 1,2
%   w_p_1, w_p_2: weight of tone model m = 1,2
% frequencies high and low in cents
% N = length(cents_stft);

hh = waitbar(0,'Estimating melody','Name','Estimating the F0...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        
F0_pdf = zeros(M,tau);

c_1_i = zeros(M,N,H);
c_2_i = zeros(M,N,H);
for h = 1:H
    c_1_i(:,:,h) = repmat(c_0_1(h),M,N);
    c_2_i(:,:,h) = repmat(c_0_2(h),M,N);
end

for t = 1:tau
    %tic
    
    S = Signal(:,t);
    S_MN = repmat(S',M,1);
    S_MNH = repmat(S_MN,[1,1,H]);
    
    waitbar(t/tau,hh)
    % re-initialise c_1
    c_1_N = c_1_i;
    c_2_N = c_2_i;
    
    if t == 1 % first frame
        % w : M*N
        w_1 = repmat(1/M/2,M,N);
        w_2 = w_1;
    end
    
    for it = 1:3 % iteration time
        
        %G_W = zeros(M,N,H); % eq 13
        
        ph_1 = c_1_N.*G_W; % eq 11
        ph_2 = c_2_N.*G_W;
        
        %for h = 1:H
        %G_W = zeros(M,N); % eq 13
        %for F = 1:M
        %xx0 = (F+Fl-1)+1200*log2(h);
        
        %G_W(F,:) = exp(-(cents_stft-xx0).^2/2/W^2)/sqrt(2*pi*W^2); % eq 13
        
        % ph_1(F,:,h) = c_1(h,F)*G_W(F,:,h); % eq 11
        % ph_2(F,:,h) = c_2(h,F)*G_W(F,:,h);
        
        %end
        %end
        
        % p M*N
        p_1 = sum(ph_1,3); % eq,10
        p_2 = sum(ph_2,3);
        
        wp_1 = w_1.*p_1;
        wp_2 = w_2.*p_2;
        
        %Norm_factor : 1*N
        Norm_factor = sum(wp_1 + wp_2); % numerator of eq 43, 44
        Norm_factor(Norm_factor == 0) = eps;
        Nf_MN = repmat(Norm_factor,M,1);
        Nf_MNH = repmat(Nf_MN,[1,1,H]);
        
        % w_ML : M*1
        w_ML_1 = sum(S_MN.*(wp_1./Nf_MN),2); % eq 43
        w_ML_2 = sum(S_MN.*(wp_2./Nf_MN),2);
        
        % w : M*N
        w_1 = repmat(w_ML_1,1,N);
        w_2 = repmat(w_ML_2,1,N);
        
        % c_ML : M*H, eq 44
        c_ML_1 = squeeze(sum(S_MNH.*(repmat(w_1,[1,1,H]).*ph_1)./Nf_MNH,2));%./repmat(w_ML_1,H,1);
        c_ML_2 = squeeze(sum(S_MNH.*(repmat(w_2,[1,1,H]).*ph_2)./Nf_MNH,2));%./repmat(w_ML_2,H,1);
        
        %c_1 = (c_ML_1)./repmat(w_ML_1,1,H);
        %c_2 = (c_ML_2)./repmat(w_ML_2,1,H);
        
        c_1 = bsxfun(@rdivide,c_ML_1',w_ML_1');
        c_2 = bsxfun(@rdivide,c_ML_2',w_ML_2');
        
        
        %normalise
        c_1 = bsxfun(@rdivide,c_1,sum(c_1))'; % faster than repmat
        c_2 = bsxfun(@rdivide,c_2,sum(c_2))';
        %c_1 = c_1./repmat(sum(c_1,2),1,H);
        %c_2 = c_2./repmat(sum(c_2,2),1,H);
        
        % 0/0 = nan, nan is due to 0/(sum = 0)
        c_1(isnan(c_1)) = 0;
        c_2(isnan(c_2)) = 0;
        
        for h = 1:H
            c_1_N(:,:,h) = repmat(c_1(:,h),1,N); % M*N*H
            c_2_N(:,:,h) = repmat(c_2(:,h),1,N);
        end
        % for debug
        %if any(any(isnan(w_1))) || any(any(isnan(w_2))) || any(any(isnan(c_1))) || any(any(isnan(c_2)))
        %    break
        %end
    end
    %toc
    %java.lang.Runtime.getRuntime.freeMemory

    
    F0_pdf(:,t) = w_1(:,1)+w_2(:,1); % eq 20
end

delete(hh)

end

