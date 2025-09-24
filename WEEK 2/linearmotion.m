function [position] = linearmotion(time)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    arguments (Input)
        time
        
    end

    arguments (Output)
        position
    end

    position = 0 + 1*time + randn; 
end