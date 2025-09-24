function [position] = parabolicmotion(time)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    time
end

arguments (Output)
    position
end

    position = 0 + 1*(time)^2  + randn;
end