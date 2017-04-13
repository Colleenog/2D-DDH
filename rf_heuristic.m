function [roc_X, roc_Y, AUC, roc_X_hold, roc_Y_hold, AUC_hold] = rf_heuristic (index_cv)

% rng(1)
% index_cv = ceil(randperm(9));
 
load adequate1;
index_ML = [75:78,89:115,126:142];
adequate_new= adequate_ML(:, index_ML);
index_ML2 = [1,2,4:7,9,11,13,14,16:21,25,26,28,32:34,36,38:43,45:47];
adequate_new= adequate_new(:, index_ML2);

load inadequate_coronal1;
inadequate_coronal_new = inadequate_coronal_ML(:, index_ML);
index_ML2 = [1,2,4:7,9,11,13,14,16:21,25,26,28,32:34,36,38:43,45:47];
inadequate_coronal_new= inadequate_coronal_new(:, index_ML2);

load inadequate_transverse1;
inadequate_transverse_new = inadequate_transverse_ML(:, index_ML);
index_ML2 = [1,2,4:7,9,11,13,14,16:21,25,26,28,32:34,36,38:43,45:47];
inadequate_transverse_new= inadequate_transverse_new(:, index_ML2);



%% learn

ind = index_cv(1:6) * floor(size(adequate_new,1)/9) -  floor(size(adequate_new,1)/9)+1;
ind_s = floor(size(adequate_new,1)/9);% size between the index
index1 = [ind(1):ind(1)+ind_s, ind(2):ind(2)+ind_s, ind(3):ind(3)+ind_s, ind(4):ind(4)+ind_s, ind(5):ind(5)+ind_s, ind(6):ind(6)+ind_s];

adequate_learn = adequate_new(index1,:);
Y = ones(size(adequate_learn,1),1);


clear ind
ind = index_cv(1:6) * floor(size(inadequate_coronal_new,1)/9) -  floor(size(inadequate_coronal_new,1)/9)+1;
ind_s = floor(size(inadequate_coronal_new,1)/9);% size between the index
index2 = [ind(1):ind(1)+ind_s, ind(2):ind(2)+ind_s, ind(3):ind(3)+ind_s, ind(4):ind(4)+ind_s, ind(5):ind(5)+ind_s, ind(6):ind(6)+ind_s];

inadequate_coronal_learn = inadequate_coronal_new(index2,:);
Y = [Y;zeros(size(inadequate_coronal_learn,1),1)];


clear ind
ind = index_cv(1:6) * floor(size(inadequate_transverse_new,1)/9) -  floor(size(inadequate_transverse_new,1)/9)+1;
ind_s = floor(size(inadequate_transverse_new,1)/9);% size between the index
index3 = [ind(1):ind(1)+ind_s, ind(2):ind(2)+ind_s, ind(3):ind(3)+ind_s, ind(4):ind(4)+ind_s, ind(5):ind(5)+ind_s, ind(6):ind(6)+ind_s];

inadequate_transverse_learn = inadequate_transverse_new(index3,:);

Y = [Y;zeros(size(inadequate_transverse_learn,1),1)];

X = [adequate_learn;inadequate_coronal_learn;inadequate_transverse_learn];


%rng(1); % For reproducibility, rng controls generation of random numbers
BaggedEnsemble = TreeBagger(50,X,Y,'OOBPred','On','OOBVarImp','On');



%% predict
clear ind
ind = index_cv(7:9) * floor(size(adequate_new,1)/9) -  floor(size(adequate_new,1)/9)+1;
ind_s = floor(size(adequate_new,1)/9);% size between the index
index1 = [ind(1):ind(1)+ind_s, ind(2):ind(2)+ind_s, ind(3):ind(3)+ind_s];

adequate_predict = adequate_new(index1,:);
GT = ones(size(adequate_predict,1),1);


clear ind
ind = index_cv(7:9) * floor(size(inadequate_coronal_new,1)/9) -  floor(size(inadequate_coronal_new,1)/9)+1;
ind_s = floor(size(inadequate_coronal_new,1)/9);% size between the index
index2 = [ind(1):ind(1)+ind_s, ind(2):ind(2)+ind_s, ind(3):ind(3)+ind_s];

inadequate_coronal_predict = inadequate_coronal_new(index2,:);
GT = [GT;zeros(size(inadequate_coronal_predict,1),1)];


clear ind
ind = index_cv(7:9) * floor(size(inadequate_transverse_new,1)/9) -  floor(size(inadequate_transverse_new,1)/9)+1;
ind_s = floor(size(inadequate_transverse_new,1)/9);% size between the index
index3 = [ind(1):ind(1)+ind_s, ind(2):ind(2)+ind_s, ind(3):ind(3)+ind_s];

inadequate_transverse_predict = inadequate_transverse_new(index3,:);
GT = [GT;zeros(size(inadequate_transverse_predict,1),1)];

X1 = [adequate_predict;inadequate_coronal_predict;inadequate_transverse_predict];



[label score cost] = predict(BaggedEnsemble,X1);

score_total = score;
label_total = label;

[roc_X, roc_Y,T,AUC] = perfcurve(GT,score_total(:,2),'1');
roc_X = interp1(roc_X,1:100); roc_Y = interp1(roc_Y,1:100);
% figure,plot(roc_X, roc_Y)
% xlabel('False positive rate'); ylabel('True positive rate')
% title('ROC for classification by logistic regression')



%% predict - holdout

% correct the feature index

load adequate1_holdout;
index_ML = [75:78,89:115,126:142];
adequate_holdout= adequate_holdout_ML(:, index_ML);
index_ML2 = [1,2,4:7,9,11,13,14,16:21,25,26,28,32:34,36,38:43,45:47];
adequate_holdout= adequate_holdout(:, index_ML2);

load inadequate_coronal1_holdout;
inadequate_coronal_holdout = inadequate_coronal_holdout_ML(:, index_ML);
index_ML2 = [1,2,4:7,9,11,13,14,16:21,25,26,28,32:34,36,38:43,45:47];
inadequate_coronal_holdout= inadequate_coronal_holdout(:, index_ML2);

load inadequate_transverse1_holdout;
inadequate_transverse_holdout = inadequate_transverse_holdout_ML(:, index_ML);
index_ML2 = [1,2,4:7,9,11,13,14,16:21,25,26,28,32:34,36,38:43,45:47];
inadequate_transverse_holdout= inadequate_transverse_holdout(:, index_ML2);

%

adequate_predict = adequate_holdout;
GT = ones(size(adequate_predict,1),1);

inadequate_coronal_predict = inadequate_coronal_holdout;
GT = [GT;zeros(size(inadequate_coronal_predict,1),1)];

inadequate_transverse_predict = inadequate_transverse_holdout;
GT = [GT;zeros(size(inadequate_transverse_predict,1),1)];

X1 = [adequate_predict;inadequate_coronal_predict;inadequate_transverse_predict];



[label score cost] = predict(BaggedEnsemble,X1);

score_total = score;
label_total = label;

[roc_X_hold, roc_Y_hold,T,AUC_hold] = perfcurve(GT,score_total(:,2),'1');
roc_X_hold = interp1(roc_X_hold,1:100); roc_Y_hold = interp1(roc_Y_hold,1:100);
% figure, plot(roc_X_hold,roc_Y_hold)
% xlabel('False positive rate'); ylabel('True positive rate')
% title('ROC for classification by logistic regression')
