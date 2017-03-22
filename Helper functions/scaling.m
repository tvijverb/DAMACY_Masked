function [ S_auto ] = scaling( S_mean, Labels, auto_mode )
%SCALING Summary of this function goes here
%   Detailed explanation goes here
S_auto = S_mean;

if auto_mode == 1
    std_whole = (std(vertcat(S_mean.Data)));
    for l1 = 1:length(S_mean)
        S_auto(l1).Data = S_mean(l1).Data./(repmat(std_whole, size(S_mean(l1).Data,1),1)); 
    end
elseif auto_mode == 2
    std_control = std(vertcat(S_mean(Labels == 0).Data));
    for l1 = 1:length(S_mean)
        S_auto(l1).Data = S_mean(l1).Data./(repmat(std_control, size(S_mean(l1).Data,1),1)); 
    end
elseif auto_mode == 3 
    for l1 = 1:length(S)
        std_indiv = std(S(l1).Data);
        S_auto(l1).Data = S_mean(l1).Data./(repmat(std_indiv ,size(S_mean(l1).Data,1),1));
    end
end

end

