function [] = SPMshade(musc, MetSPM, ymax)


% shade a figure with effect region from SPM

%% identify regions of main effect
ind = find(contains([MetSPM.colheaders], musc), 1);
crit = ind+1;
log = MetSPM.data(:,ind) > MetSPM.data(:,crit);


%% post hoc analyses
TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];

% 1 vs 3
% A = MetSPM.data(:,ind+3);
A_log = abs(MetSPM.data(:,ind+3)) > MetSPM.data(:,ind+4);

% 2 vs 3
% B = MetSPM.data(:,ind+5);
B_log = abs(MetSPM.data(:,ind+5)) > MetSPM.data(:,ind+6);

% 4 vs 3
% C = MetSPM.data(:,ind+7);
C_log = abs(MetSPM.data(:,ind+7)) > MetSPM.data(:,ind+8);

% 5 vs 3
% D = MetSPM.data(:,ind+9);
D_log = abs(MetSPM.data(:,ind+9)) > MetSPM.data(:,ind+10);


%% Loop shading and post hoc
rng = [0.3 0.2 0.1 0]; 
fnt = 22; 
width = 0.05; 
for j = 1:length(log)
    if log(j) == 1
        area([j-0.5, j+0.5], [ymax ymax], 'FaceColor',rgb('Gray'), 'EdgeColor', 'none', 'FaceAlpha',0.3);
        
        if A_log(j) == 1
                        text(j, ymax - rng(1),'.', 'Color',TotalColors(1, :), 'FontSize',fnt, 'HorizontalAlignment','center');
%             area([j-1, j+1], [ymax-rng(1)-width ymax-rng(1)+width], ...
%                 'FaceColor',TotalColors(1,:),'EdgeColor', 'none');
        end
        if B_log(j) == 1
            text(j, ymax - rng(2),'.', 'Color',TotalColors(2, :), 'FontSize',fnt, 'HorizontalAlignment','center');
%             area([j-1, j+1], [ymax-rng(2)-width ymax-rng(2)+width], ...
%                 'FaceColor',TotalColors(2,:),'EdgeColor', 'none');
        end
        if C_log(j) == 1
            text(j, ymax - rng(3),'.', 'Color',TotalColors(4, :), 'FontSize',fnt, 'HorizontalAlignment','center');
%             area([j-1, j+1], [ymax-rng(3)-width ymax-rng(3)+width], ...
%                 'FaceColor',TotalColors(3,:),'EdgeColor', 'none');
        end
        if D_log(j) == 1
            text(j, ymax - rng(4),'.', 'Color',TotalColors(5, :), 'FontSize',fnt, 'HorizontalAlignment','center');
%             area([j-1, j+1], [ymax-rng(4)-width ymax-rng(4)+width], ...
%                 'FaceColor',TotalColors(4,:),'EdgeColor', 'none');
        end
    end
end


end
