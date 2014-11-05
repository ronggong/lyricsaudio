function [ R_s ] = movaverage( hl,R )
%MOVAVERAGE moving average filter
%   construt a cardinal b-spline function (window) to low pass signal
%   hl : half length of the window
%   R  : input signal
w = [linspace(0,1,hl) linspace(1,0,hl-1)];
% padding
% R_all_t_pad = [repmat(0,hl-1,1);R_all_t;repmat(0,hl-1,1)];
 R_pad = [repmat(R(1),hl-1,1);R;repmat(R(end),hl-1,1)];
R_s = conv(R_pad,w,'valid')/sum(w);

end

