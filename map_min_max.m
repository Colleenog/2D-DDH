function [map_min, map_max] = map_min_max (probability, map)

[rows, columns] = size (map);
% getting values of map at maximum probabilities
outline = zeros(rows, columns);
for y = 1:columns-1
    [ t i] = max(probability(:,y));
     if t ~= 0
            outline(i,y) = 1;
     end        
end
    
outline_map = outline .* map;
outline_map_reshaped = reshape(outline_map, [rows*columns, 1]);

%finding mean and standard deviations
[r1,c1,map_values] = find(outline_map_reshaped);

mean_map = mean(map_values);
mean_map = min(mean_map, 0.65); mean_map = max(mean_map, 0.35); % >.35 and <.7
std_map = std(map_values); 
std_map = min(std_map, 0.25); std_map = max(std_map, 0.15);% >.15 and <.25

% finding map_min and map_max
mult_std = 1;
map_min = mean_map - 2* mult_std*std_map; %define mult_std in previous line
map_max = mean_map + .1* mult_std*std_map;

end