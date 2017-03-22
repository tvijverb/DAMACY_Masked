function [ edges ] = binedges(X, binsize, minmaxcutoff, edgesfactor)
%%
% Finds the edges for bin plots
%
% Input
%   X               = m by n data matrix with m being the cells and n the
%                     principle components
%
%   binsize1/2      = number of bins, most likely square
%   minmaxcutoff    = fraction of scores that can fall outside of biplot in
%                     most unfavourablecase for 1 measurement. default = '0.01'
%   edgesfactor     = value >1. the edges will be based on max at cutoff *
%                       edgesfactor
% Ouput
%   edges           = matrix with m by n dimensions with m being the number
%                     of PCs and n being the maximum binsize. Note the
%                     smaller bins have finite values for 1:binsize and
%                     binsize + 1: maximum binsize filled with infite
%                     values.
%%
No_PCs = size(X,2);

X_sorted = sort(X);
number_of_cells = size(X,1);
min_PC = X_sorted(max(1,ceil(0.5*minmaxcutoff*number_of_cells)),:);
max_PC = X_sorted(max(1,ceil((1-0.5*minmaxcutoff)*number_of_cells)),:);

edges = inf(No_PCs, max(binsize));
for L1 = 1:No_PCs
    edges(L1,1:binsize(L1)) = linspace( (min_PC(L1)- (edgesfactor-1)*abs(min_PC(L1)) ), ( max_PC(L1) + ( (edgesfactor-1)*abs(max_PC(L1))) ), binsize(L1) );
end
end

