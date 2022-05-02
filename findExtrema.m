function [localExtrema] = findExtrema(measurements, windowSize, funcType)
% funcType = {"max", "min"} depending on extremum desired
% returns all local extrema in "measurments" of a desired type as
% determined by "windowSize"

    buffer = zeros(100, 2);
    extremaCount = 0;
    numIter = length(measurements) - (2*windowSize - 1) + 1;
    if strcmp(funcType, "min")
        for i=1:numIter
            fail = 0;
            for j=0:(2*windowSize - 1) - 1
                if(measurements(i+windowSize) > measurements(i+j))
                    fail = 1;
                    break;
                end
            end
            if fail
                continue;
            end
            extremaCount = extremaCount + 1;    
            buffer(extremaCount, 1) = measurements(i+windowSize);
            buffer(extremaCount, 2) = i+windowSize;

        end
        
    else  % max
        for i=1:numIter
            for j=0:(2*windowSize - 1) - 1
                fail = 0;
                if(measurements(i+windowSize) < measurements(i+j))
                    fail = 1;
                    break;
                end
            end
            if fail
                continue;
            end
            extremaCount = extremaCount + 1;    
            buffer(extremaCount, :) = measurements(i+windowSize);
            buffer(extremaCount, 2) = i+windowSize;
        end
        
    end
    
    localExtrema = buffer(1:extremaCount, :);
end

