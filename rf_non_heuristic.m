function [roc_X_hog, roc_Y_hog, AUC_hog, roc_X_hold_hog, roc_Y_hold_hog, AUC_hold_hog] = rf_heuristic (index_cv)

% rng(1)
% index_cv = ceil(randperm(9));
 
load hog_adequate1;
adequate_new= hog_adequate_ML;

load hog_inadequate_coronal1;
inadequate_coronal_new = hog_inadequate_coronal_ML;

load hog_inadequate_transverse1;
inadequate_transverse_new = hog_inadequate_transverse_ML;



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
BaggedEnsemble = TreeBagger(50,X,Y,'OOBPred','On');



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

[roc_X_hog, roc_Y_hog,T,AUC_hog] = perfcurve(GT,score_total(:,2),'1');
roc_X_hog = interp1(roc_X_hog,1:100); roc_Y_hog = interp1(roc_Y_hog,1:100);
% figure,plot(roc_X_hog, roc_Y_hog)
% xlabel('False positive rate'); ylabel('True positive rate')
% title('ROC for classification by logistic regression')


%% predict - holdout

% correct the feature index

load hog_adequate1_holdout;
adequate_holdout= hog_adequate_holdout_ML;

load hog_inadequate_coronal1_holdout;
inadequate_coronal_holdout = hog_inadequate_coronal_holdout_ML;


load hog_inadequate_transverse1_holdout;
inadequate_transverse_holdout = hog_inadequate_transverse_holdout_ML;



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

[roc_X_hold_hog, roc_Y_hold_hog,T,AUC_hold_hog] = perfcurve(GT,score_total(:,2),'1');
roc_X_hold_hog = interp1(roc_X_hold_hog,1:100); roc_Y_hold_hog = interp1(roc_Y_hold_hog,1:100);
% figure, plot(roc_X_hold_hog,roc_Y_hold_hog)
% xlabel('False positive rate'); ylabel('True positive rate')
% title('ROC for classification by logistic regression')


