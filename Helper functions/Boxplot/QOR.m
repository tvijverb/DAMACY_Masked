data6 = vertcat(S6.Data)';
N_Var = length(S(1).Data(1,:));
N_ID = length(S6);
paired = 0;

for i = 1:N_Var
    y(i,:) = quantile(data6(i,:),[0.025 0.25 0.50 0.75 0.975]);
end
S_temp = S6;
S_temp2 =S6;
S_O = S6;

for j = 1:N_ID
    for i = 1:N_Var
        S_temp(j).Data(:,i)=S(j).Data(:,i) <= y(i,1);
        S_temp2(j).Data(:,i)=S(j).Data(:,i) >= y(i,5);
    end
    S_temp(j).Data2(:,1) = ~((sum(S_temp(j).Data,2)+sum(S_temp2(j).Data,2)) > 3);
    S_O(j).Data = S(j).Data(S_temp(j).Data2,:);
end