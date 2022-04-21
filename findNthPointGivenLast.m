function [prevDex, newPt] = findNthPointGivenLast(lastPt, prevDex, lenIncrement, inputArr)
    remainingLen = lenIncrement;
    while (remainingLen> norm(lastPt - inputArr(prevDex+1, :)))
        remainingLen = remainingLen - norm(lastPt - inputArr(prevDex+1, :));
        lastPt = inputArr(prevDex+1, :);
        prevDex = prevDex + 1;
    end
    
    ratio = remainingLen/norm(inputArr(prevDex+1, :) - inputArr(prevDex, :));
    newPt = lastPt + ratio*(inputArr(prevDex+1, :) - lastPt);
end