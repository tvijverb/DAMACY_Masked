function plot_2d_histograms(Histogram, edges, LDS, VAR, VariableNames, pc_used, l1, label)
%%
% Plots a histogram and PCA loadings on top and gives a title and x and y label. 
% Input:
% Histogram:        The histogram you want to plot
% edges:            A 2xbinsize matrix with the binedges
% LDS:              Loadings of PCA
% VAR:              Variance of PCA
% VariableNames:    The VariableNames of the Loadings
% pc_used:          The PCs used in the PCA model
%
% Output:
% No output but a figure will appear
%
% For example:
% plot_2d_histograms(squeeze(Histogram(x,:,:)), edges, LDS, VAR, VariableNames(var_used))
%
%
% Written by G.H. Tinnevelt on 20-2-2015 at the Radboud University Nijmegen
% Edited by G.H. Tinnevelt on 22-4-2015, able to use different PCs
%%
figure
imagesc(edges(1,:),edges(2,:),Histogram')
colormap(flipud(gray));
hold on
    title(['Histogram ID: ' num2str(l1) ' Label: ' num2str(label)], 'fontsize', 22)
    set(gca, 'YDir', 'normal');
    set(gca, 'Ylim', [-15 15], 'Xlim', [-15 15]);
    xlabel(['PC1 (' num2str(VAR(1)) '%)'])
    ylabel(['PC2 (' num2str(VAR(2)) '%)'])

DA_FLplotloadings(LDS(:,pc_used), VariableNames, [1 0 0])
set(gca, 'YDir', 'normal');
% title({'Histogram of class ' num2str(Labels(l1))}, 'fontsize', 22)
%xlabel(['PC' num2str(pc_used(1)) ' (' num2str(VAR(pc_used(1))) '%)'], 'fontsize', 18)
%ylabel(['PC' num2str(pc_used(2)) ' (' num2str(VAR(pc_used(2))) '%)'], 'fontsize', 18)
