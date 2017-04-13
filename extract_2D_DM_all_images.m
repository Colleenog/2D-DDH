clear all
close hidden all

load('trained_classifier_heuristic_features')

path = 'C:\Users\niamul\Documents\MATLAB\ddh2D\adequate_numbered';
filename = 'a (5)';
path_total = strcat(path,'\',filename, '.png');
img = imread(path_total);

for i=258:-1:1
    close all
    filename = strcat('a (', num2str(i), ').png')
    %img=imread('a (12).png');
    img=imread(filename);

    img = double (img (:,:,1) );

    % load classifier
    
    
    [im, feature, alpha, beta, fhc, I,L,A,edge_r, edge_c,center_final, radius_final]=adequacyFeature_alpha_beta_fhc(img);
    adequacy = adequacy_finder(feature,BaggedEnsemble);
    
    %display lines and circle
    x1 = edge_c; y1=edge_r;
    m_A = tand(90-A);
    c_A = -m_A * x1 + y1;
    m_L = tand(90-L);
    c_L = -m_L * x1 + y1;
    m_I = tand(90-I);
    c_I = -m_I * x1 + y1;
    straight_line = [m_A, c_A; m_I, c_I; m_L, c_L];
    show_straight_line (straight_line, im, edge_c)
    viscircles(center_final,radius_final) 
    snapnow
    
    % adequacy and dysplasia metrics

    adequacy
    alpha
    beta
    fhc

end



%%
% holdout dataset
for i=50:-1:1
    close all
    filename = strcat('holdout (', num2str(i), ').png')
    %img=imread('a (12).png');
    img=imread(filename);

    img = double (img (:,:,1) );

    % load classifier
    load('trained_classifier_heuristic_features')

    [im, feature, alpha, beta, fhc, I,L,A,edge_r, edge_c,center_final, radius_final]=adequacyFeature_alpha_beta_fhc(img);
    adequacy = adequacy_finder(feature,BaggedEnsemble);
    
    
    %display lines and circle
    x1 = edge_c; y1=edge_r;
    m_A = tand(90-A);
    c_A = -m_A * x1 + y1;
    m_L = tand(90-L);
    c_L = -m_L * x1 + y1;
    m_I = tand(90-I);
    c_I = -m_I * x1 + y1;
    straight_line = [m_A, c_A; m_I, c_I; m_L, c_L];
    show_straight_line (straight_line, im, edge_c)
    viscircles(center_final,radius_final) 
    snapnow
    
    % adequacy and dysplasia metrics

    adequacy
    alpha
    beta
    fhc

end



 %%
% % inadequate dataset
% adequacy_total = 0;
% for i=135:-1:1
%     close all
%     filename = strcat('z (', num2str(i), ').png')
%     %img=imread('a (12).png');
%     img=imread(filename);
% 
%     img = double (img (:,:,1) );
% 
%     % load classifier
%     load('trained_classifier_heuristic_features')
% 
%     [im, feature, alpha, beta, fhc, I,L,A,edge_r, edge_c,center_final, radius_final]=adequacyFeature_alpha_beta_fhc(img);
% 
%     %display lines and circle
%    
%     % adequacy and dysplasia metrics
%     adequacy = adequacy_finder(feature,BaggedEnsemble)
%     adequacy_total(i)=adequacy;
% 
% 
% end

