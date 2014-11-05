function [ chroma_fil ] = chromafil( freq_stft, lowmidi,highmidi )
%CHROMAFIL generate a chorma matrix filter
%   input: freq_stft - frequency indice of fft
%          low midi  - lowest center frequency of chroma filter
%          high mid  - highest center frequency of chroma filter
%   output: chroma_fil - matrix 12 * nfreq
midiTF = 12*log2(freq_stft / 440)+69;
nFreq = length(freq_stft);
chroma_fil = zeros(12,nFreq);

for ii = 1:12
    for n = lowmidi:highmidi
        if mod(n,12) + 1 == ii
            x = abs(n - midiTF);
            H = 1/2*tanh(pi*(1-2*x))+1/2;           % filter low passs
            chroma_fil(ii,:) = chroma_fil(ii,:) + H;
        end
    end
    
    chroma_fil(ii,:) = chroma_fil(ii,:) / sum(chroma_fil(ii,:));
end

end

