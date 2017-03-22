function [S_scores2, LDS2] = flip_svd(S_scores2,LDS2,i)
%FLIP_SVD Summary of this function goes here
%   Detailed explanation goes here

for q = 1:length(S_scores2)
        S_scores2(q).Data(:,i) = S_scores2(q).Data(:,i) ./(repmat(-1, size(S_scores2(q).Data,1),1));
end

for q = 1:length(LDS2(1,:))
    for j = 1:length(LDS2(:,1))
        LDS2(:,i) = LDS2(:,i) .* repmat(-1,[length(LDS2(:,1)) 1]);
    end
end


end

