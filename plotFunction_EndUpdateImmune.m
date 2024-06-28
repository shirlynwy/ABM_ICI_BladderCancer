function M = plotFunction_EndUpdateImmune(M,in,mean_nonzero_pd1) % storing order below is the ONLY reason I need to output M

%% immune probabilities
if ~isempty(in.p_matrix_imm)
    p_temp = in.p_matrix_imm;
    p_temp(p_temp==0) = NaN;
    [~,order] = sort(std(p_temp./mean(p_temp,1,'omitnan'),[],1,'omitnan')); %
    xs = linspace(0,1,size(in.p_matrix_imm,1));
    ys = sortrows(in.p_matrix_imm(:,order));
    for i = 1:size(p_temp,2)
        set(M.fig.imm_prob_bar(i),'XData',xs,'YData',ys(:,i),'CData',M.fig.imm_prob_colors(order(i),:))
    end

    if length(xs)>1
        M.fig.ax(M.fig.imm_prob_ind).XLim = [0,1] + .5*[-1,1]*xs(2);
    end
    M.fig.ax(M.fig.imm_prob_ind).XTick = [];
    if ~isequal(order,M.fig.imm_prob_order)
        legend(M.fig.ax(M.fig.imm_prob_ind),...
            M.fig.imm_prob_bar(size(p_temp,2):-1:1),...
            M.fig.imm_events(flip(order)),...
            'AutoUpdate','off','Location','northwest')
        M.fig.imm_prob_order = order; % storing this information is the ONLY reason I need to output M
    end
end

%% receptor plot
M.fig.receptor_plots(2).XData(end+1) = M.t;
M.fig.receptor_plots(2).YData(end+1) = mean_nonzero_pd1;