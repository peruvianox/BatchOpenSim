function S = xml_clean(S)

fNames = fieldnames(S);
nField = length(fNames);

for f = 1:nField
    if isstruct(S.(fNames{f}))
        nStruct = length(S.(fNames{f}));
        if gt(nStruct,1)
            for s=1:nStruct
                subS(s) = xml_clean(S.(fNames{f})(s));
            end
            S.(fNames{f}) = subS;
        else
            S.(fNames{f}) = xml_clean(S.(fNames{f}));
        end
    elseif ~isstruct(S.(fNames{f}))
        if strcmp(fNames{f},'COMMENT')
            S=rmfield(S,'COMMENT');
        elseif isnumeric(S.(fNames{f}))
            S.(fNames{f})=num2str(S.(fNames{f}));
        end
    end
end

