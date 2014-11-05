function [ G_W ] = G_W_calcul( pitch_int,M,cents_stft,N,H,W )
%G_W_CALCUL calculate G eq.13 Goto 2004
% M : length of search interval
% N : length of signal
% H : number of harmonics
% W : std of gaussien

xx0_right = zeros(M,N,H);
for h = 1:H
    xx0_right(:,:,h) = repmat(h,M,N);
end
xx0 = repmat(pitch_int',[1,N,H]) + 1200*log2(xx0_right);
G_W = exp(-(repmat(cents_stft,[M,1,H])-xx0).^2/2/W^2)/sqrt(2*pi*W^2); % eq 13

end

