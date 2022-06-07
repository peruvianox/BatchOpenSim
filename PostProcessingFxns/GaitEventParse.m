function[L_Times, R_Times] = GaitEventParse(L_Strikes, L_Offs, R_Strikes, R_Offs)

% calculate gait event timing as % gait cycle

%% LEFT Gait Events
Strikes = L_Strikes; 
Offs = L_Offs;
Opp_Strikes = R_Strikes; 
Opp_Offs = R_Offs; 

% trim opposite side to within foot strikes
ind1 =  Opp_Strikes > min(Strikes);
ind2 = Opp_Strikes < max(Strikes);
opp_strikes = Opp_Strikes(ind1+ind2 == 2);

% trim side to contain the strikes adjacent to opposite
strikes = Strikes; 
if opp_strikes(1) > strikes(1) && opp_strikes(1) < strikes(2)
    
else 
    strikes(1) = []; 
end

ToDel = zeros(1,length(strikes)); 
for i = 1:length(strikes)-1
    if strikes(i) < opp_strikes(i)  
        
    else
        if strikes(i+1) < opp_strikes(i)  
            ToDel(i) = 1;
        end
    end
        
    if strikes(i+1) > opp_strikes(i)
        % if strikes(end) > opp_strikes(end) && strikes(end-1) < opp_strikes(end)
        
    else
        ToDel(i+1:length(strikes)) = 1;
        break
    end
end
strikes(ToDel==1) = []; 

% trim foot off events to be within foot strikes
ind1 =  Offs > min(strikes);
ind2 = Offs < max(strikes);
offs = Offs(ind1+ind2 == 2);

ind1 =  Opp_Offs > min(strikes);
ind2 = Opp_Offs < max(strikes);
opp_offs = Opp_Offs(ind1+ind2 == 2);

    
% calculate metrics for each stride
for i = 1:length(strikes) - 1
    
    % stride time
    L_Times(i).StrideTime = strikes(i+1) - strikes(i);
    
    % opp foot off
    L_Times(i).OppFootOff = opp_offs(i) - strikes(i);
    L_Times(i).OppFootOff_p = (opp_offs(i) - strikes(i)) / L_Times(i).StrideTime;
    
    % opp foot on
    L_Times(i).OppFootOn = opp_strikes(i) - strikes(i);
    L_Times(i).OppFootOn_p = (opp_strikes(i) - strikes(i)) / L_Times(i).StrideTime;
    
    % foot off
    L_Times(i).FootOff = offs(i) - strikes(i);
    L_Times(i).FootOff_p = (offs(i) - strikes(i)) / L_Times(i).StrideTime;
    
end

clearvars Strikes Offs Opp_Strikes Opp_Offs strikes offs opp_strikes opp_offs

%% RIGHT Gait Events
Strikes = R_Strikes; 
Offs = R_Offs;
Opp_Strikes = L_Strikes; 
Opp_Offs = L_Offs; 

% trim opposite side to within foot strikes
ind1 =  Opp_Strikes > min(Strikes);
ind2 = Opp_Strikes < max(Strikes);
opp_strikes = Opp_Strikes(ind1+ind2 == 2);

% trim side to contain the strikes adjacent to opposite
strikes = Strikes; 
if opp_strikes(1) > strikes(1) && opp_strikes(1) < strikes(2)
    
else 
    strikes(1) = []; 
end

ToDel = zeros(1,length(strikes)); 
for i = 1:length(strikes)-1
    if strikes(i) < opp_strikes(i)  
        
    else
        if strikes(i+1) < opp_strikes(i)  
            ToDel(i) = 1;
        end
    end
        
    if strikes(i+1) > opp_strikes(i)
        % if strikes(end) > opp_strikes(end) && strikes(end-1) < opp_strikes(end)
        
    else
        ToDel(i+1:length(strikes)) = 1;
        break
    end
end
strikes(ToDel==1) = []; 

% trim foot off events to be within foot strikes
ind1 =  Offs > min(strikes);
ind2 = Offs < max(strikes);
offs = Offs(ind1+ind2 == 2);

ind1 =  Opp_Offs > min(strikes);
ind2 = Opp_Offs < max(strikes);
opp_offs = Opp_Offs(ind1+ind2 == 2);

    
% calculate metrics for each stride
for i = 1:length(strikes) - 1
    
    % stride time
    R_Times(i).StrideTime = strikes(i+1) - strikes(i);
    
    % opp foot off
    R_Times(i).OppFootOff = opp_offs(i) - strikes(i);
    R_Times(i).OppFootOff_p = (opp_offs(i) - strikes(i)) / R_Times(i).StrideTime;
    
    % opp foot on
    R_Times(i).OppFootOn = opp_strikes(i) - strikes(i);
    R_Times(i).OppFootOn_p = (opp_strikes(i) - strikes(i)) / R_Times(i).StrideTime;
    
    % foot off
    R_Times(i).FootOff = offs(i) - strikes(i);
    R_Times(i).FootOff_p = (offs(i) - strikes(i)) / R_Times(i).StrideTime;
    
end

%% Quality control 

Ind.ST1 = [R_Times.StrideTime] > 1.5;
Ind.ST2 = [R_Times.StrideTime] < 0;
Ind.OppOff1 = [R_Times.OppFootOff] > 10;
Ind.OppOff2 = [R_Times.OppFootOff] < 0;
Ind.OppOn1 = [R_Times.OppFootOn] > 10;
Ind.OppOn2 = [R_Times.OppFootOn] < 0;
Ind.FootOff1 = [R_Times.FootOff] > 1;
Ind.FootOff2 = [R_Times.FootOff] < 0;
Ind.ToDel = sum([Ind.ST1; Ind.ST2; Ind.OppOff1; Ind.OppOff2; ...
    Ind.OppOn1; Ind.OppOn2; Ind.FootOff1; Ind.FootOff2]) > 0; 
R_Times(Ind.ToDel) = [];


Ind.ST1 = [L_Times.StrideTime] > 1.5;
Ind.ST2 = [L_Times.StrideTime] < 0;
Ind.OppOff1 = [L_Times.OppFootOff] > 10;
Ind.OppOff2 = [L_Times.OppFootOff] < 0;
Ind.OppOn1 = [L_Times.OppFootOn] > 10;
Ind.OppOn2 = [L_Times.OppFootOn] < 0;
Ind.FootOff1 = [L_Times.FootOff] > 1;
Ind.FootOff2 = [L_Times.FootOff] < 0;
Ind.ToDel = sum([Ind.ST1; Ind.ST2; Ind.OppOff1; Ind.OppOff2; ...
    Ind.OppOn1; Ind.OppOn2; Ind.FootOff1; Ind.FootOff2]) > 0; 
L_Times(Ind.ToDel) = [];


%%

% %% Stride Time
% % duration of stride time in s
% 
% % LEFT
% for n = 1:length(L_Strikes)-1
%     L_Times(n).StrideTime = L_Strikes(n+1) - L_Strikes(n);
% end
% 
% % RIGHT
% for n = 1:length(R_Strikes)-1
%     R_Times(n).StrideTime = R_Strikes(n+1) - R_Strikes(n);
% end
% 
% 
% %% FOOT OFF (% gait cycle)
% % Find the average percent of gait cycle when foot off occurs
% 
% % Find only foot off times that are in a valid gait cycle (foot strike
% % before and after the foot off event)
% 
% % LEFT
% for n = 1:length(L_Strikes)-1
%     if L_Strikes(n) < L_Offs(n) && L_Offs(n) < L_Strikes(n+1)
%         L_Times(n).FootOff = 100*(L_Offs(n)-L_Strikes(n)) / (L_Strikes(n+1) - L_Strikes(n));
%     end
% end
% 
% % RIGHT
% for n = 1:length(R_Strikes)-1
%     if R_Strikes(n) < R_Offs(n) && R_Offs(n) < R_Strikes(n+1)
%         R_Times(n).FootOff = 100*(R_Offs(n)-R_Strikes(n)) / (R_Strikes(n+1) - R_Strikes(n));
%     end
% end
% 
% 
% % for t = 1:length(RGaitCycles(1,1,:))
% %     % Find only foot off times that are in a valid gait cycle
% %     for n = 1:length(R_Offs)
% %         if (RGaitCycles(1,1,t) < R_Offs(1,n)) && (R_Offs(1,n) < RGaitCycles(1,2,t))
% %             RFOtimes(1,n) = R_Offs(1,n);
% %         end
% %     end
% %     % Calculate percent of gait cycle when foot off occurs for each gait cycle
% %     for m = 1:length(RFOtimes)
% %         %         percent_Rfoot_offs(:,t) = (RFOtimes(1,m)-RGaitCycles(1,1,t)) / (RGaitCycles(1,2,t)-RGaitCycles(1,1,t));
% %         adj_RGaitCycles(t).FootOff  = 100*(RFOtimes(1,m)-RGaitCycles(1,1,t)) / (RGaitCycles(1,2,t)-RGaitCycles(1,1,t));
% %         adj_RGaitCycles(t).FootOffTime = RFOtimes(m);
% %     end
% % end
% 
% %% OPPOSITE FOOT OFF (% gait cycle)
% % INITIAL DOUBLE SUPPORT
% % Find the average percent of gait cycle when opposite foot off occurs
% 
% % LEFT
% for n = 1:length(L_Strikes) - 1
%     if R_Offs(1) > L_Strikes(n) && R_Offs(1) < L_Strikes(n+1)
%         Off_Ind = n;
%         break
%     end
% end
% new_L_Strikes = L_Strikes(Off_Ind:end);
% 
% if R_Offs(1) < L_Strikes(1)
%     % if right foot off occurs before first left foot contact
%     % ignore first right foot off
%     new_R_Offs = R_Offs(2:end);
% else
%     % else keep first right foot off time
%     new_R_Offs = R_Offs;
% end
% if new_R_Offs(end) > L_Strikes(end)
%     % if last right foot off times occurs after last left foot strike
%     % remove last right foot off
%     new_R_Offs = new_R_Offs(1:end-1);
% end
% GCs = min([length(new_L_Strikes), length(new_R_Offs)]);
% % if length(new_R_Offs) == 1
% %     %     percent_L_OppFootOff = ((new_R_Offs(1)-L_Strikes(1))/(L_Strikes(2)-L_Strikes(1)));
% %     L_Times(1).OppFootOff = 100*((new_R_Offs(1)-L_Strikes(1))/(L_Strikes(2)-L_Strikes(1)));
% % elseif length(new_R_Offs) == length(L_Strikes)
% %     for k = 1:length(new_R_Offs)-1
% %         %         percent_L_OppFootOff(k) = ((new_R_Offs(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
% %         L_Times(k).OppFootOff = 100*((new_R_Offs(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
% %     end
% % else
% %     for k = 1:length(new_R_Offs)
% %         %         percent_L_OppFootOff(k) = ((new_R_Offs(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
% %         L_Times(k).OppFootOff = 100*((new_R_Offs(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
% %     end
% % end
% for n = 1:GCs-1
%     L_Times(n).OppFootOff = 100*((new_R_Offs(n)-new_L_Strikes(n))/(new_L_Strikes(n+1)-new_L_Strikes(n)));
% end
% 
% % RIGHT
% for n = 1:length(R_Strikes) - 1
%     if L_Offs(1) > R_Strikes(n) && L_Offs(1) < R_Strikes(n+1)
%         Off_Ind = n;
%         break
%     end
% end
% new_R_Strikes = R_Strikes(Off_Ind:end);
% 
% if L_Offs(1) < R_Strikes(1)
%     % if left foot off occurs before first right foot contact
%     % ignore first left foot off
%     new_L_Offs = L_Offs(2:end);
% else
%     % else keep first left foot off time
%     new_L_Offs = L_Offs;
% end
% if new_L_Offs(end) > R_Strikes(end)
%     % if last left foot off occurs after last right foot strike
%     % remove last left foot off
%     new_L_Offs = new_L_Offs(1:end-1);
% end
% GCs = min([length(new_R_Strikes), length(new_L_Offs)]);
% % if length(new_L_Offs) == 1
% %     %     percent_R_OppFootOff = ((new_L_Offs(1)-R_Strikes(1))/(R_Strikes(2)-R_Strikes(1)));
% %      R_Times(1).OppFootOff = 100*((new_L_Offs(1)-R_Strikes(1))/(R_Strikes(2)-R_Strikes(1)));
% % elseif length(new_L_Offs) == length(R_Strikes)
% %     for k = 1:length(new_L_Offs)-1
% %         %         percent_R_OppFootOff(k) = ((new_L_Offs(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k)));
% %          R_Times(k).OppFootOff = 100*((new_L_Offs(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k)));
% %     end
% % else
% %     for k = 1:length(new_L_Offs)
% %         %         percent_R_OppFootOff(k) = ((new_L_Offs(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k)));
% %          R_Times(k).OppFootOff = 100* ((new_L_Offs(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k)));
% %     end
% % end
% for n = 1:GCs-1
%     R_Times(n).OppFootOff = 100*((new_L_Offs(n)-new_R_Strikes(n))/(new_R_Strikes(n+1)-new_R_Strikes(n)));
% end
% 
% %% OPPOSITE FOOT CONTACT (% gait cycle)
% % FINAL DOUBLE SUPPORT
% % Find the average percent of gait cycle when opposite foot contact occurs
% 
% % LEFT
% if R_Strikes(1) < L_Strikes(1)
%     % if first right foot contact occurs before first left foot contact
%     % ignore first right foot contact
%     new_R_Strikes = R_Strikes(2:end);
% else
%     % else keep first right foot contact time
%     new_R_Strikes = R_Strikes;
% end
% if new_R_Strikes(end) > L_Strikes(end)
%     % if last right foot strike time occurs after last left foot strike
%     % remove last right foot strike
%     new_R_Strikes = new_R_Strikes(1:end-1);
% end
% if length(new_R_Strikes) == length(L_Strikes)
%     for k = 1:length(new_R_Strikes)-1
%         %         percent_L_OppFootOn(k) = ((new_R_Strikes(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
%         L_Times(k).OppFootOn = 100*((new_R_Strikes(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
%     end
% else
%     for k = 1:length(new_R_Strikes)
%         %         percent_L_OppFootOn(k) = ((new_R_Strikes(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
%         L_Times(k).OppFootOn = 100*((new_R_Strikes(k)-L_Strikes(k))/(L_Strikes(k+1)-L_Strikes(k)));
%     end
% end
% 
% % RIGHT
% if L_Strikes(1) < R_Strikes(1)
%     % if first left foot contact occurs before first right foot contact
%     % ignore first left foot contact
%     new_L_Strikes = L_Strikes(2:end);
% else
%     % else keep first left foot strike time
%     new_L_Strikes = L_Strikes;
% end
% if new_L_Strikes(end) > R_Strikes(end)
%     % if last left foot strike time occurs after last right foot strike
%     % remove last left foot strike
%     new_L_Strikes = new_L_Strikes(1:end-1);
% end
% if length(new_L_Strikes) == length(R_Strikes)
%     for k = 1:length(new_L_Strikes)-1
%         %         percent_R_OppFootOn(k) = ((new_L_Strikes(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k))); % frames
%         R_Times(k).OppFootOn = 100*((new_L_Strikes(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k)));
%     end
% else
%     for k = 1:length(new_L_Strikes)
%         %         percent_R_OppFootOn(k) = ((new_L_Strikes(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k))); % frames
%         R_Times(k).OppFootOn = 100*((new_L_Strikes(k)-R_Strikes(k))/(R_Strikes(k+1)-R_Strikes(k)));
%     end
% end

end
