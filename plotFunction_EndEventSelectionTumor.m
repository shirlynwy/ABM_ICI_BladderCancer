function plotFunction_EndEventSelectionTumor(M, prob_matrix)

if M.NT==0 % then no tumor cells, so just delete the data
    for i = 1:2
        set(M.fig.tum_prob_bar(i),'XData',[],'YData',[])
    end
    return
end
%% tumor probabilities
order = [2,1]; % first, separate by fast vs slow, then engaged/non-engaged
xs = linspace(0,1,size(prob_matrix,1));
ys = [sortrows(prob_matrix(~M.ha_log,order));sortrows(prob_matrix(M.ha_log,order))];
for i = 1:2
    set(M.fig.tum_prob_bar(i),'XData',xs,'YData',ys(:,i))
end
if length(xs)>1
    M.fig.ax(M.fig.tum_prob_ind).XLim = [0,1] + .5*[-1,1]*xs(2);
end

if any(~M.ha_log) && any(M.ha_log)
    xts_L = xs(sum(~M.ha_log))*[0,.5,1] + xs(2)*[-.5,0,.5];
    xts_R = xs(2)/2 + [(xs(sum(~M.ha_log))+1)/2,1];

    M.fig.ax(M.fig.tum_prob_ind).XTick = [xts_L,xts_R];
    M.fig.ax(M.fig.tum_prob_ind).XTickLabel = {'|','LA','|','HA','|'};
    M.fig.ax(M.fig.tum_prob_ind).XAxis.TickLength(:) = 0;
else
    M.fig.ax(M.fig.tum_prob_ind).XTick = (xs(1)+xs(end))/2;
    if ~any(~M.ha_log)
        M.fig.ax(M.fig.tum_prob_ind).XTickLabel = 'HA';
    else
        M.fig.ax(M.fig.tum_prob_ind).XTickLabel = 'LA';
    end
end
