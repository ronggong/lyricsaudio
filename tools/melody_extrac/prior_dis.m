function [ c_0_1,c_0_2] = prior_dis%( cents_stft )
%PRIOR_DIS define the prior distribution of tone model
%   
H = 16;
h0 = 1;
U = 5.5;
%W = 17;
G = zeros(1,H);
%F = 4000; 
%N = length(cents_stft);

%%%%%% equation (46) %%%%%%%%%
c_0_1 = zeros(1,H);
c_0_2 = zeros(1,H);

for h = 1:H
    G(h) = exp(-(h-h0)^2/2/U^2)/sqrt(2*pi*U^2);
    c_0_1(h) = G(h);
end

for h = 1:2:H
    c_0_2(h) = G(h);
end

for h = 2:2:H
    c_0_2(h) = G(h)*2/3;
end

c_0_1 = c_0_1/sum(c_0_1);
c_0_2 = c_0_2/sum(c_0_2);

%%%%%%%%%%%%%%%%%%%%%%
% p_0_1 = zeros(1,N);
% p_0_2 = zeros(1,N);
% 
% for h = 1:H
%     G_W = zeros(1,N);
%     xx0 = F+1200*log2(h);
%     for ii = 1:N
%         xx = cents_stft(ii);
%         G_W(ii) = exp(-(xx-xx0)^2/2/W^2)/sqrt(2*pi*W^2);
%     end
%     
%     p_0_1 = p_0_1 + c_0_1(h)*G_W;
%     p_0_2 = p_0_2 + c_0_2(h)*G_W;   
% end

end

