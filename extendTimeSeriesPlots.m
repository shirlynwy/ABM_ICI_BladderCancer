function extendTimeSeriesPlots(M,ei)

if M.plot_pars.plot_fig
    AXES = findall(allchild(M.fig.handle),'type','Axes');
    for ax_num = 1:length(AXES)
        set(M.fig.handle,'currentaxes',AXES(ax_num))
        H = gca;
        if isempty(findall(allchild(H),'type','Line'))
            continue; % no lines on these axes so probably not plotting over time
        end

        yR = max([H.YLim(1)*(1+eps()); arrayfun(@(L) max(L.YData,[],'all'), findall(allchild(H),'type','Line'))  ]  );
        yL = max(0,H.YLim(1));
        line(M.t*[1 1],[yL yR],'DisplayName','Update Line')
        if ei<size(M.events.times,1)
            H.XLim(2) = max(H.XLim(2),M.events.times(ei+1,1));
        end
    end
end