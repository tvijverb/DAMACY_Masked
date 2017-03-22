function [Hists] = NDhist(X, Edges, smooth_factor)
% to BIN score values into equally sized bins in specified ranges
%
% Code adapted from: http://www.mathworks.nl/matlabcentral/fileexchange/13352-smoothhist2d
% Based on article 'Enhancing scatterplots with smoothed densities'. Eilers et al. (2004)
%
% Input
%   MAT            = Data, Samples x Variables. (number of bins/dims is set by
%                     ABSedges.
%   ABSedges       = Ndims x Nbins, from min to max for each dimension.
%   smooth_factor  = lambda (article), 1-5 normal, 0 for no smoothing.
%
% Ouput
%   binMAT     = [binsize x binsize x ... ] for N dimensions.
%                Histogram is normalised by the number of samples.
%                Each column of input MAT fills dimensions sequencially.
% DA - August 2014
% GT - 29 Septembre 2014: Able to accept different binsizes

% % Binning
[Hists] = sub_DA_bin(X, Edges);
% % Smoothing   (See article, Eiler et al 2004.)
dims = size(Hists); %
ndims = length(dims); %
% Function: prod([]) , reshape([],[dim end],prod[dim1-end-1])
if smooth_factor > 0;
   for L1 = 1:ndims; % Main loop, For each dimension I will unfold.
      order = 1:ndims; order(1) = L1; order(L1) = 1; % 1. for permutation of dimensions
      Hists = permute(Hists, order); % 2. bring L1 to first dimention
      Hists = Hists(:, 1:end); % 3. unfold. (don't use Reshape because it keeps last dimention rather than first (because of sub2ind))
      Hists = sub_DA_smooth1D(Hists, smooth_factor); % 4. smoothing, see subfunction
      Hists = reshape(Hists, dims(order)); % 5. reshape to permuted matrix
      Hists = permute(Hists, order); % 6. permute to original dimensions
      
      % reduced readability but less lines: 3+4+5+6 in one line.
      % binMAT=permute(reshape(sub_DA_smooth1D(binMAT(:,1:end),smooth_factor),dims(order)),order);
   end
end

end

%% ----------------------------------------------------------------------------
function [Hists] = sub_DA_bin(X, Edges)
% will create a multidimensional histogram.
BinArray = zeros(size(X, 1), size(Edges, 1)); % preallocate, each sample gets Ndim coordinates
for L1 = 1:size(Edges, 1)
    ind = isfinite(Edges(L1,:));
    ind_num = 1:size(Edges,2);
    ind_num = ind_num(ind);
    Sizes(L1) = length(ind_num);
    Edges(L1, ind) = [- Inf Edges(L1, ind_num(2:end-1)), Inf];
    [~, BinArray(:, L1)] = histc(X(:, L1), Edges(L1, :)); % Later permute dimensions 1 and 2 for imagesc, X as PC1
end

for L1 = 1:size(Edges,1)
    BinArray((BinArray(:,L1) == 2),:) = [];
    X((BinArray(:,L1) == 2),:) = [];
    BinArray(BinArray(:,L1) == Sizes(L1)-1,:) = [];
    X(BinArray(:,L1) == Sizes(L1)-1,:) = [];
end
%Sizes = repmat(size(Edges, 2), 1, size(Edges, 1));
if length(Sizes) == 1
    Sizes = [Sizes 1];
end
Hists = (accumarray(BinArray, 1, Sizes) ./ (size(X, 1))); % normalised ND histogram

end

%% -----------------------------------------------------------------------------
function Z = sub_DA_smooth1D(Y, lambda) % See source FileExchange / Article
if lambda == 0;
   Z = Y;
else
    if size(Y,1) == 1
        Y = Y';
    end
    [m, ~] = size(Y); % smoothing in 1 direction.
   E = eye(m);
   D1 = diff(E, 1); % first derivative
   D2 = diff(D1, 1); % second derivative
   P = lambda .^ 2 .* D2'*D2 + 2.*lambda .* D1' * D1; %
   Z = (E + P) \ Y; % E+P=normalized
end
end