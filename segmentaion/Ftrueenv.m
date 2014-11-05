% True-Envelope estimator
%
% DESCRIPTION
%
%  The original interface of mtltruenv.mex:
% [env,coef,envphas,freqmel,freqlin]=mtltrueenv(ampspec,order,winsize
%                                     ,mode,maxit,prec,presmooth,freqbreak)
%  Please follow this interface for any extension of this code.
%
% REFERENCES
%
% [1] Roebel, A., and Rodet, X.: Efficient Spectral Envelope Estimation And Its
%     Application To Pitch Shifting And Envelope Preservation, Proc. Digital
%     Audio Effects (DAFx), DAFX, 30-35, 2005
%
% [2] Roebel, A., and Rodet, X.: Real Time Signal Transposition With Envelope
%     Preservation In The Phase Vocoder, ICMC, 2005
%
% COPYRIGHT
%  Coded at University of Crete, Greece
%  In order to make the SVLN distribution of Ircam independent on mex functions
%
% AUTHOR
%  Gilles Degottex <degottex@csd.uoc.gr>
%
% $Id: Ftrueenv.m,v 1.1 2011/06/21 16:38:37 degottex Exp $

function [E ci, n] = Ftrueenv(S, order, winsize, mode, maxit, prec, presmooth_factor)

    %if ~isempty(winsize); disp('Use of winsize not yet implemented'); end
    if nargin<4; mode = 1; end
    if nargin<5; maxit = 10000; end
    if nargin<6; prec = 0.01; end
    dblim = log(Finvlp(prec));
    if nargin<7; presmooth_factor = 2; end

    if bitand(mode,1) % Use smoothing window
        order = round(1.2*order); % [1] 1.66, [matmtl] 1.2
        win = hamming(2*order+1);
        win = win((end-1)/2+1:end);
    end

    dftlen = length(S);

    A = log(abs(S));
    ci= zeros(1+order,1);

    if bitand(mode,512)
        pE = Ftrueenv(S, round((order)/presmooth_factor), [], 1);
        slim = round(0.25*dftlen/order);
        A(1:slim) = log(abs(pE(1:slim)));
        A(end-slim+2:end) = log(abs(pE(end-slim+2:end)));
    end
    A0 = A;

    n=0;
    max_diff = Inf;
    while n<200 && max_diff>dblim

        cc = ifft(A);

        cip = cc;
        cip = [cip(1); 2*cip(2:1+order)];

        if bitand(mode,1)               % use smoothing window
            % Eo = sqrt(sum((2*cc(order+2:end/2)).^2));
            Eo = sqrt(sum((2*cc(order+2:round(end/2))).^2));
            Ei = sqrt(sum(cip(1:order+1).^2));
            lambda = (Ei+Eo)/Ei;
            ci = lambda.*win.*(cip-ci) + ci; % [1] eq (5)
        else
            ci = cip;
        end

        lV = fft(ci, dftlen);
        A = max(A,real(lV));            % max of log amplitudes

        max_diff = max(A0-real(lV));    % can create over-shot

        n = n+1;

        if 0
            V3clear();
            V3spec(exp(A0), 1, 'k');
            V3spec(exp(A), 1, 'g');
            V3spec(exp(lV), 1, 'b'); subplot(312);xlim([0 0.04]);
        end
    end

    E = exp(lV);
    
return

% [E cc o] = Fspec_tenv(S, sr, f0, length(S), 56);
