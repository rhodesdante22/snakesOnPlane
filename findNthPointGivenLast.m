function [prevDex, newPt] = findNthPointGivenLast(lastPt, prevDex, lenIncrement, inputArr)
    remainingLen = lenIncrement;
    wrapAround = 0;
    skip = 0;
    if prevDex>=size(inputArr,1)
        wrapAround = 1;
        skip = 1;
    end
    if ~skip
        while (remainingLen> norm(lastPt - inputArr(prevDex+1, :)))
            remainingLen = remainingLen - norm(lastPt - inputArr(prevDex+1, :));
            lastPt = inputArr(prevDex+1, :);
            prevDex = prevDex + 1;
            if prevDex>=size(inputArr,1)
                wrapAround = 1;
                break
            end
        end
    end
    if wrapAround
        ratio = remainingLen/norm(lastPt - inputArr(1, :));
        if isnan(ratio)
            ratio = 0;
        end
        newPt = lastPt + ratio*(inputArr(1, :) - lastPt);
        return
    end

    if isnan(ratio)
        ratio = 0;
    end
    ratio = remainingLen/norm(inputArr(prevDex+1, :) - inputArr(prevDex, :));
    newPt = lastPt + ratio*(inputArr(prevDex+1, :) - lastPt);
end