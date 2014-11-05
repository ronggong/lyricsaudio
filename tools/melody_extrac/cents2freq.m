function [freq ] = cents2freq( cents )
%CENTS2FREQ Summary of this function goes here
freq = 2.^(cents/1200)*(440*2^(0.25-5));


end

