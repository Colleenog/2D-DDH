function adequacy = adequacy_finder(feature,BaggedEnsemble)

%% adequacy

index_ML = [75:78,89:115,126:142];
adequate_new= feature(:, index_ML);
index_ML2 = [1,2,4:7,9,11,13,14,16:21,25,26,28,32:34,36,38:43,45:47];
adequate_new= adequate_new(:, index_ML2);

[label score cost] = predict(BaggedEnsemble,adequate_new);


    adequacy=score(1,2);

  % computation here %
