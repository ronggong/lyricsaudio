function [ ssm_chro2,L ] = chromanr( ssm_chro,hopsize,fs )
%CHROMANR chroma noise remover
%   ssm_chro : chroma lag ssm
tau = length(ssm_chro);
ssm_chro2 = zeros(tau);
L = round(1.2/hopsize*fs);
% testmat = zeros(tau);
h = waitbar(0,'1','Name','ssm noise remove...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
for ii = (2*L):(tau-L+1)
    if getappdata(h,'canceling')
        break
    end
    if (mod(ii,20) == 0)
        waitbar(ii/tau,h,sprintf('%d,%d',ii,tau))
    end
    for jj = L:(ii-L+1)
        lm_l = sum(ssm_chro(jj,ii-L+1:ii))/L;       % local mean left
        lm_r = sum(ssm_chro(jj,ii:ii+L-1))/L;       % local mean right
        lm_u = sum(ssm_chro(jj:jj+L-1,ii))/L;              % local mean up
        lm_d = sum(ssm_chro(jj-L+1:jj,ii))/L;              % local mean down
        lm_ur = sum(ssm_chro(jj:jj+L-1,ii:ii+L-1))/L;              % local mean upright
        lm_dl = sum(ssm_chro(jj-L+1:jj,ii-L+1:ii))/L;              % local mean downleft
        vec = [lm_l lm_r lm_u lm_d lm_ur lm_dl];
        max_v = max(vec);
        min_v = min(vec);
        
        if (max_v == lm_l || max_v == lm_r)
            ssm_chro2(jj,ii) = ssm_chro(jj,ii) - min_v;
        else
            ssm_chro2(jj,ii) = ssm_chro(jj,ii) - max_v;
        end
        
%         ii
%         jj
%         ii-L+1
%         jj+L-1
%         jj-L+1
%         ii+L-1
%         testmat(jj,ii-L+1:ii) = 1;
%         testmat(jj,ii:ii+L-1) = 1;
%         testmat(jj:jj+L-1,ii) = 1;
%         testmat(jj-L+1:jj,ii) = 1;
%         testmat(jj:jj+L-1,ii:ii+L-1) = 1;
%         testmat(jj-L+1:jj,ii-L+1:ii) = 1;
    end
end
delete(h)

end

