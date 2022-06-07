GCData = table();
i = 1;
clc;
for s = 1:length(Subjects)
    for t = 1:length(Subjects(s).Trials)
        if contains(Subjects(s).Trials(t).name, {'static','Static'})
            continue
        end
        
        if isfield( Subjects(s).Trials(t), 'GC') == 0
            continue
        elseif isempty(Subjects(s).Trials(t).GC)
            continue
        end
        
        SubjName = Subjects(s).name; 
        TrialName = Subjects(s).Trials(t).name;
        L_GC = Subjects(s).Trials(t).GC.L.RankTimes(1:10, :);
        R_GC = Subjects(s).Trials(t).GC.R.RankTimes(1:10, :);
        
        GCData(i,:) = {SubjName, TrialName, reshape(L_GC', [1,20]), reshape(R_GC', [1,20])};
        i = i+1; 
     
    end
end
    

writetable(GCData, 'GCData.csv'); 