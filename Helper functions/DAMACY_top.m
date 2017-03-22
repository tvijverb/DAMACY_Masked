function [b_acc,n_LV,w,yhat] = DAMACY_top(X,ID,Y,numFolds,paired)
%% Created by Thomas Vijverberg on 09-06-2015 at Radboud University Nijmegen
% Adapted from Gerjen Tinneveld damacy_top_cv.m
% Last edited by Thomas Vijverberg on 09-0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data matrix (MAT) 'X' with dimension ID x (total_measurements)

% Input vector (vect) ID (unique(ID) x 1) of measurements

% Input total number of latent variables to try (

%Input if you want the resulting DAMACY PLOT (BOOL) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paired DATA
if paired
    ff = unique(ID);
    acc = zeros(numFolds, length(ff));
    for l1 = 1:length(ff)
        tmp_id = ID == ff(l1);
        Xtest = X(tmp_id,:);
        Xtrain = X(~tmp_id,:);
        Y_test = Y(tmp_id,:);
        Ytrain = Y(~tmp_id,:);
        [MC_Ytrain, mean_Y] = mncn(Ytrain);
        [MC_train, mean_X] = mncn(Xtrain);
         MC_test = Xtest - repmat(mean_X, size(Xtest,1),1);

        for i_fold = 1:numFolds
            [w(:,i_fold,l1),q,P_o,W_o] = OPLS(MC_train,MC_Ytrain,i_fold);
            [~,~,y_hat(:,l1,i_fold)] = OPLSpred(MC_test,P_o,W_o,w(:,i_fold,l1),q,i_fold);
            Y_pred1 = y_hat(:,l1,i_fold) > 0 - mean_Y;
            Y_pred2 = y_hat(:,l1,i_fold) < 0 - mean_Y;
            Y_pred = Y_pred1 - Y_pred2;
            acc(i_fold,l1) = sum(Y_pred == Y_test)/length(Y_test);
        end
    end
    
%% UNPAIRED DATA
elseif ~paired
    % Randomize data in X
    idx1 = find(Y == -1);
    idx2 = find(Y == 1);
    R1 = randperm(length(idx1));
    R2 = randperm(length(idx2));
    R1 = idx1(R1);
    R2 = idx2(R2);
    ff = length(ID)/2;
    acc = zeros(numFolds,ff);
    w=zeros(length(X),numFolds,ff);
    mean_Y = mean(Y);
    for l1 = 1:ff
        % Create training and testset
        tmp_id = zeros(length(ID),1);
        %tmp_id([l1 l1+ff]) = 1;
        tmp_id(l1) = 1;
        R_test = [R1(l1); R2(l1)];
        R_train = setdiff([R1; R2], R_test);
        X_test = X(R_test,:);
        Y_test = Y(R_test,:);
        %Y_train = Y(R_train,:) - mean_Y;
        Y_train = Y(R_train,:);
        X_train = X(R_train,:);
%       [MC_Ytrain, mean_Y] = mncn(Y_train);
%       [MC_train, mean_X] = mncn(X_train);
%       MC_test = X_test - repmat(mean_X, size(X_test,1),1);

        for i_fold = 1:numFolds
            % Run OPLS and get accuracy
            [w(:,i_fold,l1),q,P_o,W_o] = OPLS(X_train,Y_train,i_fold);
            [~,~,y_hat(:,l1,i_fold)] = OPLSpred(X_test,P_o,W_o,w(:,i_fold,l1),q,i_fold);
            Y_pred1 = y_hat(:,l1,i_fold) > 0 - mean_Y;
            Y_pred2 = y_hat(:,l1,i_fold) < 0 - mean_Y;
            Y_pred = Y_pred1 - Y_pred2;
            acc(i_fold,l1) = sum(Y_pred == Y_test)/length(Y_test);   
        end
    end
end
% See which LV had the highest accuracy and send back to parent function
    [b_acc, n_LV] = max(mean(acc'));
    w = mean(squeeze(w(:,n_LV,:)),2);
    yhat = reshape(y_hat(:,:,n_LV)',[],1);
end






