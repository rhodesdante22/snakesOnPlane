function [area] = find_area(contour)

area = 0; 
for i = 1:length(contour)-1

    if i <= length(contour)/2
        line1(i,:) = [contour(i,1)-contour(i+1, 1), contour(i,2)-contour(i+1,2)]; 
        line2(i,:) = [contour(i,1)-contour(length(contour)+1-i, 1), contour(i,2)-contour(length(contour)+1-i, 2)];
        line3(i,:) = [contour(i+1,1)-contour(length(contour)+1-i, 1), contour(i+1,2)-contour(length(contour)+1-i, 2)];
    else
        line1(i,:) = [contour(i,1)-contour(i+1, 1), contour(i,2)-contour(i+1,2)]; 
        line2(i,:) = [contour(i,1)-contour(length(contour)-i, 1), contour(i,2)-contour(length(contour)-i, 2)];
        line3(i,:) = [contour(i+1,1)-contour(length(contour)-i, 1), contour(i+1,2)-contour(length(contour)-i, 2)];
    end

    a = sqrt(line1(i,1)^2+line1(i,2)^2); 
    b = sqrt(line2(i,1)^2+line2(i,2)^2);
    c = sqrt(line3(i,1)^2+line3(i,2)^2);

    s = (a+b+c)/2; 
    area = area + sqrt(s*(s-a)*(s-b)*(s-c)); 

end