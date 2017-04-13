function [ map, karamalis_bone, combined_attenuation] = bone_feature_extraction(prep_im)

im = prep_im; 
%map
alpha = 1.0; beta = 90; gamma = 0.05;
[ map ] = confMap(im, alpha, beta, gamma); 

[rows, columns] = size (map);

%karamalis_bone
L = zeros (rows, columns);
attenuation = zeros (rows, columns);
shadowing = zeros (rows, columns);
combined = zeros (rows, columns);


g = round( size(im,2)/40)   -2; %g = size of window, 1mm

%map2= (map+0.1)/ 1.1; % small value added to avoid infinity

for i=1:rows  % -16 is for matrix dimension because of window
    for j=1:columns          
        index1 = max(i-g, 1);
        index2 = max(i-1, 1);
        index3 = max(j-g, 1);
        index4 = min(j+g, columns);
        index5 = min(i+1, rows);
        index6 = min(i+g, rows);
        
        temp_window1 = im ( index1:index2, index3:index4 ); sum1 = sum(sum(temp_window1));
        temp_window2 = map ( index5:index6, index3:index4 ); sum2 = sum(sum(temp_window2)); 
        temp_window3 = map ( index1:index2, index3:index4 ); sum3 = sum(sum(temp_window3));   

% mean based: takes more time, and not really that effective
%         temp_window1 = im ( index1:index2, index3:index4 ); sum1 =mean(mean(temp_window1)); 
%         temp_window2 = map ( index5:index6, index3:index4 ); sum2 = mean(mean(temp_window2)); 
%         temp_window3 = map ( index1:index2, index3:index4 ); sum3 = mean(mean(temp_window3)); 
        
        attenuation (i,j) =  sum3-sum2;
        shadowing (i,j) =  sum3/(sum2+.1); 
        
        L (i,j) =  sum1 - sum2; 
    end
end
karamalis_bone = max(L, -0);
karamalis_bone = circshift(karamalis_bone, [-g 0]);
for i=(rows-g):(rows)
    karamalis_bone(i,:) = karamalis_bone(rows-g-1,:);
end
karamalis_bone = normalize (karamalis_bone);

attenuation = normalize(attenuation);
shadowing = normalize(shadowing);
combined_attenuation = attenuation + shadowing;
combined_attenuation = normalize(combined_attenuation);


end