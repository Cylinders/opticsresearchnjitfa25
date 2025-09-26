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
    complex = times((0.1 * linearmotion(t) + randn), exp(linearmotion(t)*1i * k));

end