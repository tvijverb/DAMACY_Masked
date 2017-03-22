%% DAMACY script
%% This script uses FLOOD and DAMACY combined
% 
% Written by G.H. Tinnevelt on 26-2-2015 at Radboud University Nijmegen
% Last edited by Thomas Vijverberg on 10-06-2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input data structure (STRUCT) 'S' with dimension 1 x (total_measurements)
    %Must contain atleast contain:
    %Data input (S.Data)
    %IDs input (S.ID)
    %Labels S.Labels
    
%Input VariableNames matrix (MAT) with dimension 1 x (number_variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables Used
% Cleaning workspace except required variables
clearvars -except S VariableNames

% Variable default settings -- description for each is presented in script
% below --
paired_data = 1;
center_mode = 3;
auto_mode = 1;
pc_used = [1 2];
N_PCs = length(pc_used);
N_Var = length(S(1).Data(1,:));
N_ID = length(S);
binsize = repmat(200,1,N_PCs); 
PLS_iterations = 5;
Labels = vertcat(S.Labels);
ID = vertcat(S.ID);
smooth_factor = 1;
w_rescaling = (1/3);
w_cutoff = 1;

%% Preprocessing
%% outlier removal
% Calculate 95% quantiles (5% considered outlier) for each ID,variable.
% If cell of individual has 4 or more variables considered outlier, it is
% removed. This is done to correct instrumental errors in FACS data.

data = vertcat(S.Data)';
for i = 1:N_Var
    y(i,:) = quantile(data(i,:),[0.025 0.25 0.50 0.75 0.975]);
end
S_temp = S;
S_temp2 =S;
S_O = S;

for j = 1:N_ID
    for i = 1:N_Var
        S_temp(j).Data(:,i)=S(j).Data(:,i) <= y(i,1);
        S_temp2(j).Data(:,i)=S(j).Data(:,i) >= y(i,5);
    end
    S_temp(j).Data2(:,1) = ~((sum(S_temp(j).Data,2)+sum(S_temp2(j).Data,2)) > 2);
    S_O(j).Data = S(j).Data(S_temp(j).Data2,:);
end

clear data
%% meancenter
fprintf('\n------------------Meancentering data------------------\n');
tic
%MEANCENTERING
%1 mean center over all controls
%2 median center over all controls
%3 mean center over all individuals 
%4 median center over all individuals

if paired_data == 1
    if center_mode == 5 || center_mode == 6
        S_mean = unpaired_centering(S_O,center_mode -2);
    elseif sum(center_mode == 1:4) == 1
        S_mean = unpaired_centering(S_O, center_mode);
    else
        S_mean = S_O;
    end
else
    if sum(center_mode == 1:4) == 1
        S_mean = unpaired_centering(S_O, center_mode);
    else
        S_mean = S_O;
    end
end 

toc
%% scaling
fprintf('\n------------------Scaling data------------------\n');
tic
% SCALING
%1 Pareto scaling complete dataset
%2 UV scaling over control group
%3 UV scaling over individuals
S_auto = scaling(S_mean, Labels, auto_mode);

clear S_mean auto_mode
toc

%% Create PCA model on control data with number of PCs given
fprintf('\n------------------Creating PCA model------------------\n');
tic
%Principle Component Analysis and selecting --pc_used-- components
% eg. pc_used = [1 2 3];

data = vertcat(S_auto.Data);
[SCR, LDS, VAR] = pcafunction(data);

for L1=1:(min(length(S(1).Data(1,:)),7))
        disp([ '| Expained Variance - PC ' num2str(L1) ' = ' num2str(VAR(L1)) ' %']);
end 

for l1 = 1:length(S_auto)
    S_scores(l1,1).Data = S_auto(l1).Data*LDS(:,pc_used);
end

clear S_auto center_mode L1 l1 data

toc
%% Plot PCAfigures
% for iis= [1 46]
%     %1:size(ID,1)
%     figure;
%   
%     plot(S_scores(iis).Data(:,1), S_scores(iis).Data(:,2), '*');
%     hold on
%     DA_FLplotloadings(LDS(:,1:2),VariableNames, [1 0 0])
%     set(gca, 'YDir', 'normal');
%     title([ 'ID: '   num2str(ID(iis)) ' Label: ' num2str(Labels(iis))], 'fontsize', 22)
%     xlabel(['PC1 (' num2str(VAR(1)) '%)'], 'fontsize', 18)
%     ylabel(['PC2 (' num2str(VAR(2)) '%)'], 'fontsize', 18)
% end

%% Create binedges
fprintf('\n------------------Create bins histogram------------------\n');
tic
%Create histogram edge grid based on the maximum and minimum values over 
%each column. 

X = vertcat(S_scores.Data);
minimum_value = min(X);
maximum_value = max(X);
for i = 1:length(minimum_value)
    edges(i,:) = linspace(minimum_value(i),maximum_value(i),binsize(1));
end

clear X
toc
%% Create Histograms
fprintf('\n------------------Filling bins histogram------------------\n');
tic
%Filling Histograms with scaled data, the smooth factor spreads the peaks
%in the histogram. 0 = no smoothing (1-5) are respectively slight smoothing
%to heavy smoothing. 
Histograms = zeros([length(S) binsize]);
for l1 = 1:N_ID
    Histograms(l1,:,:,:,:,:,:) = NDhist(S_scores(l1).Data, edges, smooth_factor); 
     %plot_2d_histograms(squeeze(Histograms(l1,:,:)), edges, LDS, VAR, VariableNames, pc_used,ID(l1),Labels(l1))
end

clear l1
toc
%% remove empty bins (DAMACY)
%Note: Empty bins should not be taken into account to circumvent seperation
%based on empty vs not empty. The empty bins are not 0, but have a very low
%value due to the smoothing.
fprintf('\n------------------Remove empty bins------------------\n');
tic

X = reshape(Histograms, size(Histograms,1), prod(binsize));
ind = std(X) > 10^-6;
X = X(:,ind);

toc
%% create class vector Y and train and test set (DAMACY)
fprintf('\n------------------Create trainingset and testset------------------\n');
tic
%Preparing data for regression
Y2 = Labels == 0;
Y1 = ~Y2;
Y = Y1 - Y2;
mean_Y = mean(Y);
X = X - repmat(mean(X), size(X,1), 1);
Ym = Y - mean(Y);

clear Y2 Y1 l1
toc
%% OPLS (DAMACY)
fprintf('\n------------------Perform OPLS-DA------------------\n');
tic
%Regression using OPLS-DA, J. Tryg, S. Wold; Orthogonal projections to
%latent structures (OPLS), Umea University Sweden

[b_acc,n_LV,w,yhat] = DAMACY_top(X,ID,Y,PLS_iterations,paired_data);

%Reshape OPLS weights back to histogram space
w_box = zeros(binsize);
w_box(ind) = w;

clear X Ym q P_o W_o Y_pred Y_pred1 Y_pred2
toc
%% Backtrace coordinates of weights in PCA space
fprintf('\n------------------Calculate coordinates of weights in PCA space------------------\n');
tic
% Find which bin each cell occupies, multiple occupants per cell
% possible

for q = 1:length(S_scores)
    S_scores(q).score_edges = zeros(length(S_scores(q).Data),N_PCs);
    for z = 1:N_PCs
        for j = 1:length(S_scores(q).Data)
            element = find( edges(z,:) >= S_scores(q).Data(j,z));
            S_scores(q).score_edges(j,z) = element(1);
        end
    end
end

toc

%% Get weights in PCA space
fprintf('\n------------------Get weights in PCA space------------------\n');
tic
% Scale weights to minimize the overexpression of healthy cells
% Get for each cell in struct (S) its corresponding OPLS weight
w_box_edited = w_box;
 w_box_edited = real(w_box.^(w_rescaling));
  stdv_wbe = std(reshape(w_box_edited,[1 ],[]));
 w_box_edited = (w_box_edited > w_cutoff*stdv_wbe);
 
 if N_PCs == 2;
     for i = 1:length(S_scores)
         S_scores(i).PCAweights = zeros(length(S_scores(i).Data),1);
         for j = 1:length(S_scores(i).Data)
             S_scores(i).PCAweights(j,1) = w_box_edited(S_scores(i).score_edges(j,1),S_scores(i).score_edges(j,2));
         end
     end
 elseif N_PCs == 3;
     for i = 1:length(S_scores)
         S_scores(i).PCAweights = zeros(length(S_scores(i).Data),1);
         for j = 1:length(S_scores(i).Data)
             S_scores(i).PCAweights(j,1) = w_box_edited(S_scores(i).score_edges(j,1),S_scores(i).score_edges(j,2),S_scores(i).score_edges(j,3));
         end
     end
 end

clear w_box_edited
toc
%% Scale weights respecting total amount of cells
fprintf('\n------------------Rescale weights------------------\n');
tic
% Respect original amount of cells, normalization

multiplier(1:length(S_scores)) = 0;
for i = 1:length(S_scores)
         multiplier(i) = abs(length(S_scores(i).PCAweights(:,1))/sum(S_scores(i).PCAweights(:,1)));
         S_scores(i).PCAweights = repmat(multiplier(i), [size(S_scores(i).PCAweights) 1]) .* S_scores(i).PCAweights;
end

toc
%% Project weights on original scores matrix in PCA space
% fprintf('\n------------------Project weights on original scores matrix------------------\n');
% tic
% % Apply weights to each row, one cell, in scores matrix
% 
% for i = 1:length(S_scores)
%     S_scores(i).PCAweighted = S_scores(i).Data .* repmat(S_scores(i).PCAweights,1,N_PCs);
% end
% 
% clear abs_w_box
% 
% toc
%% Project weights on original matrix
fprintf('\n------------------Project weights on original matrix------------------\n');
tic
% Apply weights to each row, one cell, in original matrix
S_mean_new = unpaired_centering(S_O, 1);
S_auto_new = scaling(S_mean_new, Labels, 1);
S_new = S_auto_new;

for i = 1:length(S_scores)
    S_new(i).Data = S_auto_new(i).Data .* repmat(S_scores(i).PCAweights,[1 length(S(1).Data(1,:))]);
    S_new(i).Labels = S(i).Labels;
end

clear abs_w_box S_mean_new S_auto_new
toc
%% Get new loadings
fprintf('\n------------------Calculate new loadings------------------\n');
tic
%Principle Component Analysis and selecting --pc_used-- components
% eg. pc_used = [1 2 3];
S_scores2 = S_scores;

[~,LDS2,VAR2] = pcafunction(vertcat(S_new.Data));

for l1 = 1:length(S_new)
    S_scores2(l1,1).Data = S_new(l1).Data*LDS2(:,pc_used);
    S_scores2(l1,1).Data( ~any(S_scores2(l1,1).Data,2), : ) = [];
end
%LDS = LDS(:,pc_used);

clear S_auto center_mode L1 l1
toc
%% Sign Check SVD
fprintf('\n------------------Sign check SVD------------------\n');
tic
%Check if the scores and loading have not been flipped by the SVD algorithm
%computational shortcoming of the matlab SVD algorithm. Gives both scores
%and loadings backflip if they are flipped.

for i = 1:length(pc_used)
    correct = 0;
    for j = 1:length(LDS(:,1))
        correct = correct + ((LDS2(j,i) > 0) == (LDS(j,i) > 0));
    end
    if correct <(length(LDS(:,1))/2)
        [S_scores2, LDS2] = flip_svd(S_scores2,LDS2,i);
    end 
end

toc
%% Create Histograms
fprintf('\n------------------Filling bins histogram------------------\n');
tic
%Create histogram edge grid based on the maximum and minimum values over 
%each column. 
%Filling Histograms with scaled data, the smooth factor spreads the peaks
%in the histogram. 0 = no smoothing (1-5) are respectively slight smoothing
%to heavy smoothing.

X = vertcat(S_scores2.Data);
minimum_value = min(X);
maximum_value = max(X);
for i = 1:length(minimum_value)
    edges2(i,:) = linspace(minimum_value(i),maximum_value(i),binsize(1));
end

Histograms2 = zeros([length(S) binsize]);
for l1 = 1:length(S)
    Histograms2(l1,:,:,:,:,:,:) = NDhist(S_scores2(l1).Data, edges2, smooth_factor); 
%     plot_2d_histograms(squeeze(Histograms(l1,:,:)), edges, LDS, VAR, VariableNames, pc_used)
end

clear minmaxcutoff edgesfactor X
toc
%% remove empty bins (DAMACY)
%Note: Empty bins should not be taken into account to circumvent seperation
%based on empty vs not empty. The empty bins are not 0, but have a very low
%value due to the smoothing.
fprintf('\n------------------Remove empty bins------------------\n');
tic

X = reshape(Histograms2 , size(Histograms2 ,1), prod(binsize));
ind = std(X) > 10^-6;
X = X(:,ind);

toc
%% create class vector Y and train and test set (DAMACY)
fprintf('\n------------------Create trainingset and testset------------------\n');
tic
% Preparing data for OPLS regression
Y2 = Labels == 0;
Y1 = ~Y2;
Y = Y1 - Y2;
mean_Y = mean(Y);
X = X - repmat(mean(X), size(X,1), 1);
Ym = Y - mean(Y);

clear Y2 Y1 l1
toc
%% OPLS (DAMACY)
fprintf('\n------------------Perform OPLS-DA 2------------------\n');
tic
%Regression using OPLS-DA, J. Tryg, S. Wold; Orthogonal projections to
%latent structures (OPLS), Umea University Sweden

[b_acc2,n_LV2,w2,yhat2] = DAMACY_top(X,ID,Y,PLS_iterations,paired_data);

w_box2 = zeros(binsize);
w_box2(ind) = w2;

clear Ym q2 P_o2 W_o2 Y_pred Y_pred1 Y_pred2 S_scores S_scores2 SCR
toc

%% Visualize 2 or 3 PCs weighted vector 1(DAMACY)
fprintf('\n------------------Visualise OPLS-DA result 1------------------\n');
tic
%Plotting function
if N_PCs == 2
    DAMACY_weights_2d( yhat,w_box,LDS,Y,edges,binsize,VAR,VariableNames,ID,b_acc);
elseif N_PCs == 3
    DAMACY_weights_3d( yhat,w_box,LDS,Y,edges,binsize,VAR,VariableNames,ID);
%     giffy('original weights.gif');
end
toc
%% Visualize 2 or 3 PCs weighted vector 2(DAMACY)
fprintf('\n------------------Visualise OPLS-DA result 2------------------\n');
tic
%Plotting function
if N_PCs == 2
    DAMACY_weights_2d( yhat2,w_box2,LDS2,Y,edges2,binsize,VAR2,VariableNames,ID,b_acc2);
elseif N_PCs == 3
    DAMACY_weights_3d( yhat2,w_box2,LDS2,Y,edges2,binsize,VAR2,VariableNames,ID);
%     giffy('new weights.gif');
end
toc

%% Convert scores to structure
% ID=vertcat(S.ID);
% S_scores = struct('Data', []);
% indv_sizes = zeros(length(S),1);
% 
% for i = 1: length(S)
%     S_scores(i,1).Data = scores(1:length(S(i).Data),pc_used);
%     scores(1:length(S(i).Data),:) = [];
% end

%% Plot PCAfigures
% for iis=4
%     %1:size(ID,1)
%     figure;
%   
%     plot(S_scores(iis).Data(:,1), S_scores(iis).Data(:,2), '*');
%     hold on
%     DA_FLplotloadings(LDS(:,1:2),VariableNames, [1 0 0])
%     set(gca, 'YDir', 'normal');
%     title([ 'ind '   num2str(ID(iis)) ], 'fontsize', 22)
%     xlabel(['PC1 (' num2str(VAR(1)) '%)'], 'fontsize', 18)
%     ylabel(['PC2 (' num2str(VAR(2)) '%)'], 'fontsize', 18)
%    
%     
% end

%% Applying OPLS-DA weights to histograms and plot new histograms
% Weighted_Histograms = Histograms;
% 
% for l1 = 1:length(S)
%     Weighted_Histograms(l1,:,:,:,:) = abs(squeeze(Histograms(l1,:,:,:,:)).*w_box);
% %      plot_2d_histograms(squeeze(Weighted_Histograms(l1,:,:)), edges, LDS, VAR, VariableNames, pc_used)
% end

%% Different smoothing function
% for i = 1:N_ID
%     %Histograms(i,:,:,:,:,:,:,:) = smoothn(Histograms(1,:,:,:,:,:,:,:),0.45,'robust');
%     yy = smoothn(squeeze(Histograms(i,:,:,:,:,:,:,:)),0.4,'robust');
%     Histograms(i,:,:,:,:,:,:,:) = yy;
% end
