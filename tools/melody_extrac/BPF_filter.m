function [ filter_curve,cutoff_curve_bins ] = BPF_filter( freq_stft, cutoff_curve )
%BPF_FILTER generate the filter_curve

% freq_stft: freq bins - > freq
% cutoff_curve_bins : frequency bins of 4 cut off points
% cutoff_curve: 4 numbers in Hz, bandpass_low0, bandpass_low1, high1, high0

N = length(freq_stft);
cutoff_curve_bins = zeros(1,4);
filter_curve = zeros(1,N);

for ii = 1:4
    [~,cutoff_curve_bins(ii)] = min(abs(freq_stft-cutoff_curve(ii)));
end

for ii = cutoff_curve_bins(1):cutoff_curve_bins(2)
    %filter_curve(ii) = (ii - cutoff_curve_bins(1))/(cutoff_curve_bins(2)-cutoff_curve_bins(1));
    filter_curve(ii) = sqrt(1-((ii-cutoff_curve_bins(2))/(cutoff_curve_bins(2)-cutoff_curve_bins(1)))^2);
end

filter_curve(cutoff_curve_bins(2)+1:cutoff_curve_bins(3)) = 1;

for ii = cutoff_curve_bins(3)+1:cutoff_curve_bins(4)
    filter_curve(ii) = 1 - (ii - cutoff_curve_bins(3))/(cutoff_curve_bins(4)-cutoff_curve_bins(3));
end


end

