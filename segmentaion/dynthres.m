function [ th ] = dynthres( lsf_s )
%DYNTHRES dynamic threshold eq.8
%   lsf_s : input
l_lsf = length(lsf_s);
min_start = max(mean(lsf_s),0);
    l_th = floor((max(lsf_s)-min_start)/0.01);
    th_s = zeros(l_th,1);
    for jj = 1:l_th
        th = min_start + 0.01*(jj-1);
        w1 = sum(lsf_s>=th)/l_lsf;
        w2 = 1-w1;
        u1 = mean(lsf_s(lsf_s>=th));
        u2 = mean(lsf_s(lsf_s<th));
        th_s(jj,1) = w1*w2*(u1-u2)^2;
    end
    
    [~,th] = max(th_s);
    th = min_start + 0.01*(th-1);

end

