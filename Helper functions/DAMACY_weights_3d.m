function DAMACY_weights_3d( yhat,w_box,LDS,Y,edges,binsize,VAR,VariableNames,ID)
%DAMACY_WEIGHTS_3D Summary of this function goes here
%   Detailed explanation goes here
    idx1 = Y == -1;
    idx2 = Y == 1;
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2, 3 ,[1 4] , 'Ylim', [-2 2]);hold on
    plot(0,yhat(idx1),'r*')
    text(repmat(0.25,sum(idx1),1), yhat(1:sum(idx1)), num2str(ID(1:sum(idx1))), 'VerticalAlignment','top','HorizontalAlignment','left', 'Fontsize', 5, 'color','b')
    plot(0,yhat(idx2), 'b*')
    text(repmat(0.25,sum(idx2),1), yhat(sum(idx1)+1:length(Y)), num2str(ID(sum(idx1)+1:length(Y))), 'VerticalAlignment','top','HorizontalAlignment','left', 'Fontsize', 5, 'color','r')
    line([-1 1], [0 0])
    title('yhat of controls (red) vs diseased (blue)', 'fontsize', 22)

    %w_box high is class 1
    %w_box low is class 0
    subplot(2, 3 ,[2,6]);
       load('cmap_bluered.mat')
    w_box_high = zeros(binsize);
    w_box_low = zeros(binsize);
    
    w_box_high(w_box > 0) = w_box(w_box > 0);
    w_box_low(w_box < 0) = w_box(w_box < 0);
    vol3d('cdata',permute(w_box_high, [2 1 3]), 'Xdata', edges(1,:),'Ydata', edges(2,:),'Zdata', edges(3,:),'mycolor',flipud(cmap(49:64,:)),'texture','2D');
    hold on
    vol3d('cdata',permute(abs(w_box_low), [2 1 3]), 'Xdata', edges(1,:),'Ydata', edges(2,:),'Zdata', edges(3,:),'mycolor',flipud(cmap(1:16,:)),'texture','2D');
    set(gca, 'YDir', 'normal');
    title('Loadings', 'fontsize', 22)
    campos([57.1841 -176.7572  100.6733]);
    xlabel(['PC1 (' num2str(VAR(1)) '%)'])
    ylabel(['PC2 (' num2str(VAR(2)) '%)'])
    zlabel(['PC3 (' num2str(VAR(3)) '%)'])
    DA_FLplotloadings(LDS(:,1:3),VariableNames, [1 0.5 1])

end

