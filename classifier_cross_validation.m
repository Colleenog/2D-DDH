clear all

%% heuristic
tic
for i = 1:50
    rng(i)
    index_cv = ceil(randperm(9));
    [roc_x(:,:,i), roc_y(:,:,i), auc(i), roc_x_hold(:,:,i), roc_y_hold(:,:,i), auc_hold(i)] = rf_heuristic (index_cv);
end

heuristic = mean (auc)
heuristic_hold = mean (auc_hold)


%% non heuristic - HOG
    
for i = 1:25
    rng(i)
    index_cv = ceil(randperm(9));
    [roc_x_hog(:,:,i), roc_y_hog(:,:,i), auc_hog(i), roc_x_hold_hog(:,:,i), roc_y_hold_hog(:,:,i), auc_hold_hog(i)] = rf_non_heuristic (index_cv);
end

nonheuristic = mean (auc_hog)
nonheuristic_hold = mean (auc_hold_hog)

toc
