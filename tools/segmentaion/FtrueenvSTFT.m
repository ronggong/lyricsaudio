function [ specFilter_v_abs,freq_v ] = FtrueenvSTFT( signal, sr, fftlen,hop,F0,draw_env,label )
%   FTRUEENVSTFT calculate the short-term True Envelope
%   input   signal
%   output  specFilter_v_abs
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Set default arguments %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 6
        draw_env = 0;
    end
    
    if nargin < 7
        label = '';
    end
    
    %default.sr = 48000;
    default.sr = sr;
    default.window = 'hanning';
    default.wintime = 0.03; % window size minimal: 4*T0 = 4*(1/147) = 27ms
    default.hoptime = 0.005;
    default.preemph = 0.97;
    default.midear = 0;
    default.fbtype = 'htkmel';
    default.constamp = 1;
    default.nband = 40;
    default.scale = 'lin';
    default.compress = 'log';
    default.ncep = 13;
    default.lift = 0;
    default.lwin_dmfcc = 9;
    default.lwin_ddmfcc = 5;
   
     param = default;
     
    % window size in samples
%     Lw = round(param.wintime * param.sr);
% 	Lw = min([Lw, numel(signal)]);
%   	% ensure odd length window
% 	if Lw == round(param.wintime * param.sr);
% 	Lw = Lw + ~rem(Lw,2);
% 	else
%     Lw = Lw - ~rem(Lw,2);
% 	end
%     % hop size in samples
%     hop = round(param.hoptime * param.sr);
    
    
    
    % fft size
    %nfft = 1024; % 2*2.^nextpow2(Lw)
    % extract frames of signal
    
    % only do 
    if length(signal) >= fftlen
        [signal, ~] = Fenframe(signal, fftlen , hop, 'truncate');
    end
    
    % cepstrum order
    %F0 = 147;
    ccOrder = param.sr/(2*F0);
    scaleFactor = 1.3; % 0.5 1.5
    
    % compute FFT (on the whole signal)
    %[spec_v freq_v] = Ffft(signal_v, Fs);
    
    % compute cepstrum and CC filter
    
    for ii = 1 : size(signal,2)
        
        % compute FFT (on the whole signal)
        [spec_v,freq_v] = Ffft(signal(:,ii), param.sr,fftlen);
        [specFilter_v(:,ii),cc_v(:,ii)]=Ftrueenv(spec_v, ccOrder*scaleFactor);
        
        specFilter_v(isnan(specFilter_v(:,ii)),ii) = eps;
        
        % tracer chaque enveloppe
        if draw_env ~= 0
            
            folderPath_true_env_image = '/Users/gong/Documents/MATLAB/true_env_image/';
            mkdir(folderPath_true_env_image);
            if ii == ceil (size(signal,2)/2)
                figure
                plot(freq_v, 20*log10(spec_v));
                hold on;
                plot(freq_v, 20*log10(abs(specFilter_v(:,ii))), 'r-', 'LineWidth', 1);
                % h(2)=plot(freq_v, 20*log10(abs(specFilter_v_adapt)), 'g', 'LineWidth', 1);
                % legend(h,'without adaptation','with adaptation')
                xlabel('frequency [Hz]');
                ylabel('power density spectrum [dB]');
                xlim([0 12000]);
                title(['signal et True Envelope, voyelle=' label ' f0=' num2str(F0) ' frame=' num2str(ii)])
                filename = sprintf('signal_true_env_%u_%s',ii,label);
                saveas(gca,fullfile(folderPath_true_env_image,filename),'epsc2');
                close(gcf);
            end
        end
    end
    specFilter_v_abs = abs(specFilter_v);
    
end

