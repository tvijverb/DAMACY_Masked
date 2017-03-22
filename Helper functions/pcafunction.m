function [SCR, LDS, VAR]=pcafunction(X)
%%
%Geeft de score, loading en de var matrices.
%%

[U D V]=svd(X,'econ');
SCR = U*D;
LDS = V;
VAR=diag(D.*D);
VAR = round(100*VAR/sum(VAR));
end