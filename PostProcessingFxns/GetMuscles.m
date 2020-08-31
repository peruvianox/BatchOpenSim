function [Muscles, Joints] = GetMuscles(Parsed, colheaders, type)

% organize columns of muscle outputs from opensim into left and right
% structures

%% Update Colheaders if necessary
if strcmp(type, 'met')
    if contains(colheaders{1}, '_0_')
        Str2Rep = 'metabolics_0_';
    elseif contains(colheaders{1}, '_1_')
        Str2Rep = 'metabolics_1_';
    end
    NewColHeaders = strrep(colheaders,Str2Rep, '');
    
else
    NewColHeaders = colheaders;
end
ColInds = ones(1,length(NewColHeaders));


%% Loop through column headers to assign muscles
for i = 1:length(NewColHeaders)
    
    if ColInds(i) == 0
        continue % if right muscle already found, skip left side
    end
    
    if strcmp(NewColHeaders{i}(end-1:end), '_r')
        Muscles(i).name = NewColHeaders{i}(1:end-2);
        OppInd = find(strcmp(NewColHeaders,strcat(Muscles(i).name, '_l')));
        % right muscles
        Muscles(i).data.right_col = i;
        Muscles(i).data.right_parsed = Parsed(:,i,:);
        
        % left muscles
        Muscles(i).data.left_col = OppInd;
        Muscles(i).data.left_parsed = Parsed(:,OppInd,:);
        
        ColInds(OppInd) = 0;
        
    elseif strcmp(NewColHeaders{i}, 'TOTAL')
        % save global or total metabolics
        Muscles(i).name = NewColHeaders{i};
        Muscles(i).data.global_parsed = Parsed(:,i);
        Muscles(i).data.col = i;
        
    elseif strcmp(NewColHeaders{i}, 'BASAL')
        % save global or total metabolics
        Muscles(i).name = NewColHeaders{i};
        Muscles(i).data.global_parsed = Parsed(:,i);
        Muscles(i).data.col = i;
    end
end
Joints = [];

ToDel = zeros(1, length(Muscles)); 
for i = 1:length(Muscles)
    if isempty(Muscles(i).name) % if no-name, set to delete
        ToDel(i) = 1;
    end
end

ToDel = logical(ToDel); % Delete unnecessary rows
Muscles(ToDel) = [];

%% define muscles crossing each joint
if strcmp(type, 'met')
    % split biarticular muscles?
    SplitBi = 'Yes';
    
    HipMuscles = {'glut_', 'semi','bifemlh','sar','add_','tfl','pect','grac','iliacus','psoas',...
        'quad_fem','gem','peri','rect_fem'};
    KneeMuscles = {'bifemlh','bifemsh','rect_fem','gas', 'vas_', 'semi'};
    AnkleMuscles = {'gas','soleus','tib','per_','ext_', 'flex'};
    TrunkMuscles = {'ercspn','intobl','extobl'};
    
    for i = 1:length(Muscles)
        Muscles(i).hip = 0;
        Muscles(i).knee = 0;
        Muscles(i).ankle = 0;
        Muscles(i).trunk = 0;
        
        if isempty(Muscles(i).name) % if no-name, set to delete
            ToDel(i) = 1;
        else
            ToDel(i) = 0;
            
            % dont split biarticular muscles - add to multiple joints
            if contains(Muscles(i).name, HipMuscles)
                Muscles(i).hip = 1;
            end
            if contains(Muscles(i).name, KneeMuscles)
                Muscles(i).knee = 1;
            end
            if contains(Muscles(i).name, AnkleMuscles)
                Muscles(i).ankle = 1;
            end
            if contains(Muscles(i).name, TrunkMuscles)
                Muscles(i).trunk = 1;
            end
            
            
            if strcmp(SplitBi, 'Yes')
                
                data = importdata('MuscleParameters.csv');
                
                if contains(Muscles(i).name, 'semi') || contains(Muscles(i).name, 'bifemlh')
                    Ind = find(contains(data.textdata(:,1), 'Hamstrings')) - 1;
                    h = abs(data.data(Ind, 7)); % column 7 - hip sagittal plane
                    k = abs(data.data(Ind, 10)); % column 10 - hip sagittal plane
                    t = h + k;
                    Muscles(i).hip = h / t;
                    Muscles(i).knee = k / t;
                end
                
                if contains(Muscles(i).name, 'rect_fem')
                    Ind = find(contains(data.textdata(:,1),  'Rectus femoris')) - 1;
                    h = abs(data.data(Ind, 7)); % column 7 - hip sagittal plane
                    k = abs(data.data(Ind, 10)); % column 10 - hip sagittal plane
                    t = h + k;
                    Muscles(i).hip = h / t;
                    Muscles(i).knee = k / t;
                end
                
                if contains(Muscles(i).name, 'gas')
                    Ind = find(contains(data.textdata(:,1),  'Gastrocnemius')) - 1;
                    a = abs(data.data(Ind, 11)); % column 7 - hip sagittal plane
                    k = abs(data.data(Ind, 10)); % column 10 - hip sagittal plane
                    t = a + k;
                    Muscles(i).ankle = a / t;
                    Muscles(i).knee = k / t;
                end
                
            end
            
        end
        
    end
    
    %     ToDel = logical(ToDel); % Delete unnecessary rows
    %     Muscles(ToDel) = [];
    
    
    %% Define Ankle, Knee, and Ankle variables by summing all muscles
    HipInd = 1;
    KneeInd = 1;
    AnkleInd = 1;
    TrunkInd = 1;
    TotalInd = 1;
    
    for i = 1:length(Muscles)
        if Muscles(i).hip ~= 0 % hip muscles
            Hip.Columns{HipInd,1} = Muscles(i).name;
            Hip.left_parsed(:,HipInd) = Muscles(i).data.left_parsed .* Muscles(i).hip;
            Hip.right_parsed(:,HipInd) = Muscles(i).data.right_parsed .* Muscles(i).hip;
            HipInd = HipInd + 1;
        end
        if Muscles(i).ankle ~= 0 % ankle muscles
            Ankle.Columns{AnkleInd,1} = Muscles(i).name;
            Ankle.left_parsed(:,AnkleInd) = Muscles(i).data.left_parsed .* Muscles(i).ankle;
            Ankle.right_parsed(:,AnkleInd) = Muscles(i).data.right_parsed .* Muscles(i).ankle;
            AnkleInd = AnkleInd + 1;
        end
        if Muscles(i).knee ~= 0 % knee muscles
            Knee.Columns{KneeInd,1} = Muscles(i).name;
            Knee.left_parsed(:,KneeInd) = Muscles(i).data.left_parsed .* Muscles(i).knee;
            Knee.right_parsed(:,KneeInd) = Muscles(i).data.right_parsed .* Muscles(i).knee;
            KneeInd = KneeInd + 1;
        end
        if Muscles(i).trunk ~= 0 % trunk muscles
            Trunk.Columns{TrunkInd,1} = Muscles(i).name;
            Trunk.left_parsed(:,TrunkInd) = Muscles(i).data.left_parsed .* Muscles(i).trunk;
            Trunk.right_parsed(:,TrunkInd) = Muscles(i).data.right_parsed.* Muscles(i).trunk;
            TrunkInd = TrunkInd + 1;
        end
        if strcmp(Muscles(i).name, 'TOTAL')
            Total.interp = Muscles(i).data.global_parsed;
            
%              Total.Columns{TotalInd,1} = Muscles(i).name;
%             Total.left_parsed(:,TotalInd) = Muscles(i).data.left_parsed;
%             Total.right_parsed(:,TotalInd) = Muscles(i).data.right_parsed;
%             TotalInd = TotalInd + 1;
            
        end
    end
    
    % Sum muscle energy costs across each joint
    Hip.left_total_parsed = sum(Hip.left_parsed, 2);
    Hip.right_total_parsed = sum(Hip.right_parsed, 2);
    Ankle.left_total_parsed = sum(Ankle.left_parsed, 2);
    Ankle.right_total_parsed = sum(Ankle.right_parsed, 2);
    Knee.left_total_parsed = sum(Knee.left_parsed, 2);
    Knee.right_total_parsed = sum(Knee.right_parsed, 2);
    Trunk.left_total_parsed = sum(Trunk.left_parsed, 2);
    Trunk.right_total_parsed = sum(Trunk.right_parsed, 2);
    
    Ind = contains({Muscles.name}, 'TOTAL');
    Total.parsed = Muscles(Ind).data.global_parsed;
    
    Total.left_total_parsedSum = sum([Hip.left_total_parsed,  Knee.left_total_parsed, ...
         Ankle.left_total_parsed,  Trunk.left_total_parsed], 2);
     Total.right_total_parsedSum = sum([Hip.right_total_parsed,  Knee.right_total_parsed, ...
         Ankle.right_total_parsed,  Trunk.right_total_parsed], 2);
    
    clearvars AnkleInd HipInd KneeInd TrunkInd OppInd MetStart MetEnd i
    
    %% Save joints structure
    Joints.Hip = Hip;
    Joints.Ankle = Ankle;
    Joints.Knee = Knee;
    Joints.Trunk = Trunk;
    Joints.Total = Total;
    
end

end
