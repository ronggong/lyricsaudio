function [ cents ] = freq2cents( freq )
%FREQ2CENTS Summary of this function goes here
%   Detailed explanation goes here
cents = 1200*log2(freq/(440*2^(0.25-5)));


end