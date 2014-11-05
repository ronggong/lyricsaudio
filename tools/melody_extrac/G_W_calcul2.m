function [ G_W ] = G_W_calcul2( pitch_int,M,cents_stft,N,H,W )
%G_W_CALCUL calculate G eq.13 Goto 2004
% M : length of search interval
% N : length of signal
% H : number of harmonics
% W : std of gaussien

xx0_right = zeros(N,M,H);
for h = 1:H
    xx0_right(:,:,h) = repmat(h,N,M);
end
xx0 = repmat(pitch_int,[N,1,H]) + 1200*log2(xx0_right);
G_W = exp(-(repmat(cents_stft',[1,M,H])-xx0).^2/2/W^2)/sqrt(2*pi*W^2); % eq 13

end

