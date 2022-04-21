function len = getLen(inputArr)
    len = 0;
    numPoints = size(inputArr, 1);
    for i=2:numPoints
        len = len + norm(inputArr(i, :) - inputArr(i-1, :));
    end
end