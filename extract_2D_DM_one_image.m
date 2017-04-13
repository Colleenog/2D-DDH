clear all
close hidden all

load('trained_classifier_heuristic_features')

path = 'C:\Users\niamul\Documents\MATLAB\ddh2D_femoral_head_bug_addressed\adequate_numbered\';
filename = 'a (5)';
path_total = strcat(path,'\',filename, '.png');
img = imread(path_total);
img = double (img (:,:,1) );

    % load classifier
    load('trained_classifier_heuristic_features')    
    
    % extract features and alpha, beta
    [im, feature, alpha, beta, fhc, I,L,A,edge_r, edge_c,center_final, radius_final]=adequacyFeature_alpha_beta_fhc(img);
    
    % extract adequacy
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



