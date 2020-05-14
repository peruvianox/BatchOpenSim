function[data] = ClarifyHeaders(TRCfile, data)

% get column hearders from TRC file
% headers returned back in structure

%% Clarify Headers to get acutal columns of marker positions
fid = fopen(TRCfile, 'r');
FData = textscan(fid,'%q'); 

% locate columns of interest
ColHeadStart = find(strcmp(FData{1,1},'Frame#'),1); 
ColHeadEnd = find(strcmp(FData{1,1},'X1'),1); 
i = 1; 
while isempty(ColHeadEnd) % if no X1, look for 
    i = i + 1;
    Srch = strcat('X',num2str(i));
    ColHeadEnd = find(strcmp(FData{1,1},Srch),1);
    if i == 120 % stop searching at 120 if no X# found
        ColHeadEnd = find(strcmp(FData{1,1},'1'),1); 
        break
    end
end

% move XYZ col headers down one row
ColHeaders = FData{1,1}(ColHeadStart+2: ColHeadEnd-1)'; 
if isfield(data,'colheaders')
    for i = 1:length(data.colheaders)
        data.colheaders{2,i} = data.colheaders{1,i};
    end
else
    data.colheaders{1,1} = 'Frame#';
    data.colheaders{1,2} = 'Time';
    j = 3; 
    for i = 1:length(ColHeaders)
        data.colheaders{1,j} = ColHeaders{i};
        j = j + 3;
    end
end

% copy in marker label column headers
Col = 3:3:size(data.colheaders, 2);
for i = 1:length(ColHeaders)
    data.colheaders{1,Col(i)} = ColHeaders{i};
    data.colheaders{1,Col(i) + 1} = {' '};
    data.colheaders{1,Col(i) + 2} = {' '};
end

% rename certain headers
for i = 1:length(ColHeaders)
    if strcmp(data.colheaders{1,i}, 'L.Acromium')
        data.colheaders{1,i} = 'L.Acr'; 
    end
    if strcmp(data.colheaders{1,i}, 'R.Acromium')
        data.colheaders{1,i} = 'R.Acr';
    end
end