function [ ssm_chro ] = chromassm( chrogram )
%CHROMASSM build time lag chroma similarity matrix
%   input   
%           chrogram : chroma gram
%   output  ssm_chro : similarity matrix
h = waitbar(0,'1','Name','ssm chroma calculating...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
tau = size(chrogram,2);
ssm_chro = zeros(tau);
denom = sqrt(12);
for ii = 1:tau
    if getappdata(h,'canceling')
        break
    end
    if (mod(ii,20) == 0)
        waitbar(ii/tau,h,sprintf('%d,%d',ii,tau))
    end
    for jj = 0:ii-1
        ssm_chro(jj+1,ii) = 1 - sqrt(sum((chrogram(:,ii-jj)/max(chrogram(:,ii-jj)) - chrogram(:,ii)/max(chrogram(:,ii))).^2))/denom;
    end
end
delete(h)


end

