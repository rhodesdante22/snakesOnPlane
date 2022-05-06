function longAxis = findLongAxis(points)
    maxDist = 0;
    longAxis = zeros(2);
    for i=1:size(points, 1) - 1
        for j=i+1:size(points, 1)
            if (norm(points(i, :) - points(j,:)) > maxDist)
                maxDist = norm(points(i, :) - points(j,:));
                longAxis(1,:) = points(j, :);
                longAxis(2,:)= points(i, :);
            end   
        end   
    end
end