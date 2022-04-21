function outArr = xInterp(inArr, x)
    outArr = zeros(x,2);
    len = getLen(inArr);
    lenIncrement = len/x;
    bottomLV = inArr(1,:);
    outArr(1,:) = bottomLV;
    lastDex = 1;
    for i=2:x
       [lastDex, outArr(i, :)] = findNthPointGivenLast(outArr(i-1, :), lastDex, lenIncrement, inArr);
    end    
end