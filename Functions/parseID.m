function [ParsedStride] = parseID(data, Start, End)
    % get time region
    time = data.Header; 
    startInd = find(time > Start, 1);
    endInd = find(time > End, 1);
    
    % resample
    stride = table2array(data(startInd:endInd, :)); 
    ParsedStride = resample(stride, 100, length(stride)); 
end