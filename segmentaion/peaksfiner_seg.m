function [ peaks_R_T,R_T_h_s ] = peaksfiner_seg( R_T_h, K_s_p )
%PEAKSFINER_SEG find peaks by smooth differentiel
%   R_T_h : input signal
%   K_s_p : half length of moving average window
tau = length(R_T_h);
R_T_h_s = zeros(tau,1);
 R_T_h_pad = [repmat(R_T_h(1),K_s_p,1);R_T_h;repmat(R_T_h(end),K_s_p,1)];
    for jj = 1:tau
        R_T_h_s(jj) = R_T_h_pad(jj:jj+2*K_s_p)'*(-K_s_p:K_s_p)';
        
    end
    
    peaks_R_T = [];
    for jj= 2:tau
        if R_T_h_s(jj) < 0 && R_T_h_s(jj-1) >0
            peaks_R_T(end+1) = jj-1;
        end
    end

end

