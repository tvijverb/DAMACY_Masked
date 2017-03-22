function [ S_new ] = Get_Tube( S,tube )
%GET_TUBE Summary of this function goes here
%   Detailed explanation goes here

j= 1;

for i = 1:length(S)
    if S(i).Tube == tube
        S_new(j) = S(i);
        j = j+1;
    end
end


end

