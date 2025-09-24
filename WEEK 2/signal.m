function [complex] = signal(r,t)
%SIGNAL Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    r
    t
end

arguments (Output)
    complex
end
    k = 1;
    complex = exp(r*i * k);

end