close all
clear all

%% adequate 
cd ..
cd('adequate_numbered')

for k = 1:258
    pngFileName = strcat('a ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );   
    img = imresize(img, [400, 400]);
    features = extractHOGFeatures(img);
    features = squeeze(features);
    ML(:,:,k)= double(features);
end
ML=squeeze(ML);
ML=ML';
hog_adequate_ML = ML;
cd ..
cd('classification')
save ('hog_adequate1', 'hog_adequate_ML')
cd ..

%% inadequate - coronal 
clear ML
cd('inadequate_coronal_numbered')

for k = 1:137
    pngFileName = strcat('z ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );   
    img = imresize(img, [400, 400]);
    features = extractHOGFeatures(img);
    features = squeeze(features);
    ML(:,:,k)= double(features);
end
ML=squeeze(ML);
ML=ML';
hog_inadequate_coronal_ML = ML;
cd ..
cd('classification')
save ('hog_inadequate_coronal1', 'hog_inadequate_coronal_ML')
cd ..

%% inadequate - transverse 
clear ML
cd('inadequate_transverse_numbered')

for k = 1:298
    pngFileName = strcat('z ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );   
    img = imresize(img, [400, 400]);
    features = extractHOGFeatures(img);
    features = squeeze(features);
    ML(:,:,k)= double(features);
end
ML=squeeze(ML);
ML=ML';
hog_inadequate_transverse_ML = ML;
cd ..
cd('classification')
save ('hog_inadequate_transverse1', 'hog_inadequate_transverse_ML')
cd ..

%% adequate - holdout 
clear ML
cd('adequate_numbered')

for k = 1:68
    pngFileName = strcat('holdout ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );   
    img = imresize(img, [400, 400]);
    features = extractHOGFeatures(img);
    features = squeeze(features);
    ML(:,:,k)= double(features);
end
ML=squeeze(ML);
ML=ML';
hog_adequate_holdout_ML = ML;
cd ..
cd('classification')
save ('hog_adequate1_holdout', 'hog_adequate_holdout_ML')
cd ..

%% inadequate - coronal - holdout
clear ML
cd('inadequate_coronal_numbered')

for k = 1:38
    pngFileName = strcat('holdout ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );   
    img = imresize(img, [400, 400]);
    features = extractHOGFeatures(img);
    features = squeeze(features);
    ML(:,:,k)= double(features);
end
ML=squeeze(ML);
ML=ML';
hog_inadequate_coronal_holdout_ML = ML;
cd ..
cd('classification')
save ('hog_inadequate_coronal1_holdout', 'hog_inadequate_coronal_holdout_ML')
cd ..

%% inadequate - transverse - holdout
clear ML
cd('inadequate_transverse_numbered')

for k = 1:68
    pngFileName = strcat('holdout ('...
        , num2str(k), ')','.png');
    img = imread(pngFileName);
    img = double (img (:,:,1) );   
    img = imresize(img, [400, 400]);
    features = extractHOGFeatures(img);
    features = squeeze(features);
    ML(:,:,k)= double(features);
end
ML=squeeze(ML);
ML=ML';
hog_inadequate_transverse_holdout_ML = ML;
cd ..
cd('classification')
save ('hog_inadequate_transverse1_holdout', 'hog_inadequate_transverse_holdout_ML')
cd ..