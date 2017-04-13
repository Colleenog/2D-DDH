function fhc = extract_fhc(ilium_angle, edge_of_ilium, radius, center)

x1 = edge_of_ilium(1,1);
y1 = edge_of_ilium(1,2);

x2=150;
y2= y1 + tan((ilium_angle-90).*pi/180) .* (x2-x1);

% d = distance between center and ilium line, positive if center is above
coefficients1 = polyfit([x1, x2], [y1, y2], 1); % ilium
m=(coefficients1(1,1)); C=(coefficients1(1,2));
A = m; B=-1; % Ax+By+C = 0
x1 = center(1,1); y1= center(1,2);
d = (A*x1+B*y1+C)./sqrt(A^2+B^2);

% fhc
fhc = (radius-d)/(2*radius) *100;