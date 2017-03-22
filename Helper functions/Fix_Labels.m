function [ S ] = Fix_Labels( S )
%FIX_LABELS converts strings identifying healthy IDs to binary expression

for i = 1:length(S)
    if strcmp(S(i).Labels,'"normal"')
        S(i).Labels = 0;
    else
        S(i).Labels = 1;
    end
end

