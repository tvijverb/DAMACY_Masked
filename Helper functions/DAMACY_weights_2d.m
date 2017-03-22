function DAMACY_weights_2d( yhat,w_box,LDS,Y,edges,binsize,VAR,VariableNames,ID,b_acc)
%DAMACY_WEIGHTS_2D Summary of this function goes hereDetailed explanation goes here
    idx1 = Y == -1;
    idx2 = Y == 1;
    yhat = yhat - mean([mean(yhat(idx2)) mean(yhat(idx1))]);
    figure
    set(gca,'units','normalized','outerposition',[0 0 1 1])
    subplot(2, 3 ,[1 4] , 'Ylim', [-2 2]);hold on
    plot(0,yhat(idx1),'r*')
    text(repmat(0.25,sum(idx1),1), yhat(idx1), num2str(ID(1:sum(idx1))), 'VerticalAlignment','top','HorizontalAlignment','left', 'Fontsize', 5, 'color','r')
    plot(0,yhat(idx2), 'b*')
    text(repmat(0.25,sum(idx2),1), yhat(idx2), num2str(ID(sum(idx1)+1:length(Y))), 'VerticalAlignment','top','HorizontalAlignment','left', 'Fontsize', 5, 'color','b')
    line([-1 1], [0 0])
    title('yhat of controls (red) vs diseased (blue)', 'fontsize', 22)

    subplot(2, 3 ,[2,6]);
    imagesc(edges(1,:),edges(2,:),w_box',[- 0.016 0.016]);
    hold on
    load('cmap_bluered.mat')
    set(gca, 'YDir', 'normal');
    colormap(cmap)
    format bank
    title(['OPLS weights, binsize ' num2str(binsize(1)) ', ' num2str(round(b_acc,2)) ' accuracy (LTOCV)'], 'fontsize', 22)
    xlabel(['PC1 (' num2str(VAR(1)) '%)'])
    ylabel(['PC2 (' num2str(VAR(2)) '%)'])

    hold on
    DA_FLplotloadings(LDS(:,1:2),VariableNames) %VariableNames should be first made

end

