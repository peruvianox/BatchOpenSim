function S = cleanTree(S)

fNames =  fieldnames(S);
nField = length(fNames);

for i = 1:nField
    if isstruct(S.(fNames{i}))
        nStruct = length(S.(fNames{i}));
        if gt(nStruct,1)
            for j=1:nStruct
                S.(fNames{i})(j)=cleanTree(S.(fNames{i})(j));
            end
        else
            S.(fNames{i})=cleanTree(S.(fNames{i}));
        end
    elseif strcmp(fNames{i},'COMMENT')
        S=rmfield(S,'COMMENT');
    elseif isnumeric(S.(fNames{i}))
        S.(fNames{i})=num2str(S.(fNames{i}));
    end
end


