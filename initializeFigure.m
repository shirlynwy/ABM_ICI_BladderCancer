function M = initializeFigure(M)

M.fig.handle = figure;
set(M.fig.handle,'units','normalized','outerposition',[0 0 1 1])

M.fig.nrows = 5;
M.fig.ncols = 6;

%% indices for the subplots in ax
M.fig.scatter_ind = 1;
M.fig.tum_density_ind = 2;
M.fig.imm_density_ind = 3;
M.fig.isf_density_ind = 4;
M.fig.population_ind = 5;
M.fig.receptor_ind = 6;
M.fig.drug_ind = 7;
M.fig.mass_action_ind = 8;
M.fig.event_ind = 9;
M.fig.subpop_ind = 10;
M.fig.tum_prob_ind = 11;
M.fig.imm_prob_ind = 12;
M.fig.prop_cleared_ind = 13;
M.fig.clearance_ha_ind = 14;
M.fig.clearance_la_ind = 15;
M.fig.aPD1aPDL1_ind = 16;
M.fig.killing_rate_per_tumor_ind = 17;
M.fig.cell_slice_ind = 18;
M.fig.interaction_ind = 19;
M.fig.apd1_density_ind = 20;
M.fig.afgfr3_density_ind = 21;

ind_names = fieldnames(M.fig);
ind_names = ind_names(endsWith(ind_names,'ind'));

n_axes = 0;
for i = 1:length(ind_names)
    n_axes = max(n_axes,M.fig.(ind_names{i}));
end
M.fig.ax = gobjects(n_axes,1);

%% set up axes
cellfun(@(f) evalin('caller',[f ' = M.fig.' f ';']), fieldnames(M.fig));

scatter_locs = [1,ncols+2];
cell_slice_locs = [2*ncols+1,3*ncols+2];

tumor_density_locs = 0*ncols + 3;
immune_density_locs = 1*ncols + 3;
isf_density_locs = 2*ncols + 3;
afgfr3_density_locs = 4*ncols + 1;
apd1_density_locs = 4*ncols + 2;

population_locs = 0*ncols + (5:6);
receptor_locs = 2*ncols + (5:6);
drug_locs = 3*ncols + (5:6);
% apd1apdl1_locs = 3*ncols+(3:4);
mass_action_locs = nrows*ncols;
event_locs = 3*ncols + (3:4);
interaction_locs = 4*ncols + (3:4);
subpop_locs = 1*ncols + (5:6);

tum_probs_locs = 0*ncols + 4;
imm_probs_locs = [1,2]*ncols + 4;

prop_cleared_locs = 6;
clearance_rates_la_locs = 6+ncols;
clearance_rates_ha_locs = 6+2*ncols;

killing_rate_per_tumor_locs = 5+4*ncols;

%% scatter plot
tumor_colors = winter(4);
tumor_colormap = flipud(winter);
immune_colors = autumn(2);
unengaged_color = [1 1 1];
deactivated_color = [0 0 0];
immune_colormap = flipud(autumn);
isf_colormap = copper(512);
isf_colormap(1:256,:) = [];
if M.plot_pars.plot_location
    M.fig.ax(scatter_ind) = subplot(nrows,ncols,scatter_locs);
    M.fig.ax(scatter_ind).Box = 'off';
    M.fig.ax(scatter_ind).NextPlot = 'add';
    M.fig.ax(scatter_ind).Color = M.fig.ax(scatter_ind).Parent.Color;
    M.fig.ax(scatter_ind).XTick = [];
    M.fig.ax(scatter_ind).YTick = [];
    M.fig.ax(scatter_ind).ZTick = [];
    M.fig.ax(scatter_ind).XColor = 'none';
    M.fig.ax(scatter_ind).YColor = 'none';
    M.fig.ax(scatter_ind).ZColor = 'none';

    M.fig.scatter_plots(1) = scatter3([],[],[],90,tumor_colors(1,:),'o','filled',...
        'DisplayName','Low Antigen (LA) Tumor Cells');
%     M.fig.scatter_plots(2) = scatter3([],[],[],90,tumor_colors(2,:),'o','filled',...
%         'DisplayName','LA Mutant Tumor Cells');
    M.fig.scatter_plots(3) = scatter3([],[],[],90,tumor_colors(3,:),'o','filled',...
        'DisplayName','High Antigen (HA) Tumor Cells');
%     M.fig.scatter_plots(4) = scatter3([],[],[],90,tumor_colors(4,:),'o','filled',...
%         'DisplayName','HA Mutant Tumor Cells');
    M.fig.scatter_plots(5) = scatter3([],[],[],40,immune_colors(1,:),'o','filled',...
        'DisplayName','Fast-killing CTLs');
    M.fig.scatter_plots(6) = scatter3([],[],[],40,immune_colors(2,:),'o','filled',...
        'DisplayName','Slow-killing CTLs');
    M.fig.scatter_plots(7) = scatter3([],[],[],40,unengaged_color,'o','filled',...
        'DisplayName','Unengaged CTLs','MarkerFaceAlpha',0.2);
    M.fig.scatter_plots(8) = scatter3([],[],[],40,deactivated_color,'o','filled',...
        'DisplayName','Deactivated CTLs');

    legend(M.fig.ax(M.fig.scatter_ind),'Location','SouthWest','AutoUpdate','off','Color',"none")

    M.fig.ax(M.fig.scatter_ind).Legend.Position(1) = M.fig.ax(M.fig.scatter_ind).Position(1) - M.fig.ax(M.fig.scatter_ind).Legend.Position(3);
    axis(M.fig.ax(M.fig.scatter_ind),[1,M.grid.size(1),1,M.grid.size(2),1,M.grid.size(3)] - repelem(M.grid.center,1,2))

    M.plot_pars.m = sqrt(sum(M.grid.center.^2));
    view(M.fig.ax(M.fig.scatter_ind), M.plot_pars.m*[0.8660,0,0.5])

    axis square
end

%% slice plots
warning('off','MATLAB:contour:NonFiniteData')
%% cell slice plot
M.fig.ax(cell_slice_ind) = subplot(nrows,ncols,cell_slice_locs);
M.fig.ax(cell_slice_ind).NextPlot = 'add';
M.fig.cell_slice_plot(1) = scatter3([],[],[],60,tumor_colors(1,:),'o','filled',...
    'DisplayName','Low Antigen (LA) Tumor Cells');
% M.fig.cell_slice_plot(2) = scatter3([],[],[],60,tumor_colors(2,:),'o','filled',...
%     'DisplayName','LA Mutant Tumor Cells');
M.fig.cell_slice_plot(3) = scatter3([],[],[],60,tumor_colors(3,:),'o','filled',...
    'DisplayName','High Antigen (HA) Tumor Cells');
% M.fig.cell_slice_plot(4) = scatter3([],[],[],60,tumor_colors(4,:),'o','filled',...
%     'DisplayName','HA Mutant Tumor Cells');
M.fig.cell_slice_plot(5) = scatter3([],[],[],60,immune_colors(1,:),'o','filled',...
    'DisplayName','Fast-killing CTLs','MarkerFaceAlpha',0.75);
M.fig.cell_slice_plot(6) = scatter3([],[],[],60,immune_colors(2,:),'o','filled',...
    'DisplayName','Slow-killing CTLs','MarkerFaceAlpha',0.75);
M.fig.cell_slice_plot(7) = scatter3([],[],[],60,unengaged_color,'o','filled',...
    'DisplayName','Unengaged CTLs','MarkerFaceAlpha',0.75);
M.fig.cell_slice_plot(8) = scatter3([],[],[],60,deactivated_color,'o','filled',...
    'DisplayName','Deactivated CTLs','MarkerFaceAlpha',0.75);
[~,M.fig.cell_slice_plot(9)] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
M.fig.ax(cell_slice_ind).Title.String = 'Cell Slice (z=z_{mid})';
M.fig.ax(cell_slice_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(cell_slice_ind).YLabel.String = 'Y coord (cells)';

colormap(M.fig.ax(cell_slice_ind),isf_colormap);

axis square

if ~M.plot_pars.plot_location
    legend(M.fig.cell_slice_plot(1:8),'Location','bestoutside','AutoUpdate','off')
end

%% density projection plots
%% tumor density plot
M.fig.ax(tum_density_ind) = subplot(nrows,ncols,tumor_density_locs);
[~,M.fig.tum_density_plot] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
M.fig.ax(tum_density_ind).Title.String = 'Tumor Projected onto (X,Y)';
M.fig.ax(tum_density_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(tum_density_ind).YLabel.String = 'Y coord (cells)';
M.fig.tum_density_colorbar = colorbar;

colormap(M.fig.ax(tum_density_ind),tumor_colormap);

axis square

%% immune density plot
M.fig.ax(imm_density_ind) = subplot(nrows,ncols,immune_density_locs);
[~,M.fig.imm_density_plot] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
M.fig.ax(imm_density_ind).Title.String = 'Immune Projected onto (X,Y)';
M.fig.ax(imm_density_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(imm_density_ind).YLabel.String = 'Y coord (cells)';
M.fig.imm_density_colorbar = colorbar;

colormap(M.fig.ax(imm_density_ind),immune_colormap);

axis square

%% ISF density plot
M.fig.ax(isf_density_ind) = subplot(nrows,ncols,isf_density_locs);
[~,M.fig.isf_density_plot] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
M.fig.ax(isf_density_ind).Title.String = "ISF Projected onto (X,Y)";
M.fig.ax(isf_density_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(isf_density_ind).YLabel.String = 'Y coord (cells)';
M.fig.isf_colorbar = colorbar;

colormap(M.fig.ax(isf_density_ind),isf_colormap);

axis square

%% aFGFR3 density plot
afgfr3_colormap = bone();
M.fig.ax(afgfr3_density_ind) = subplot(nrows,ncols,afgfr3_density_locs);
[~,M.fig.afgfr3_density_plot] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
M.fig.ax(afgfr3_density_ind).Title.String = "aFGFR3 Slice (z=z_{mid})";
M.fig.ax(afgfr3_density_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(afgfr3_density_ind).YLabel.String = 'Y coord (cells)';
M.fig.afgfr3_colorbar = colorbar;

colormap(M.fig.ax(afgfr3_density_ind),afgfr3_colormap);

axis square

%% aPD1 density plot
apd1_colormap = bone();
M.fig.ax(apd1_density_ind) = subplot(nrows,ncols,apd1_density_locs);
[~,M.fig.apd1_density_plot] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
M.fig.ax(apd1_density_ind).Title.String = "aPD1 Slice (z=z_{mid})";
M.fig.ax(apd1_density_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(apd1_density_ind).YLabel.String = 'Y coord (cells)';
M.fig.apd1_colorbar = colorbar;

colormap(M.fig.ax(apd1_density_ind),apd1_colormap);

axis square

warning('on','MATLAB:contour:NonFiniteData')

%% time series plots
[xtick_vals,order] = sort(M.events.times);
symbols = {'I','^','.','C'};
vals = [0,1,2,Inf];
xx = [0 max(eps(),xtick_vals(1))];

[symbol_ind,~] = find((vals==M.events.event_index(order))');
xtick_labels = symbols(symbol_ind);

for i = length(xtick_vals):-1:2
    if xtick_vals(i)==xtick_vals(i-1)
        xtick_labels{i-1} = ';';
        xtick_labels(i) = [];
    end
end
xtick_vals = unique(xtick_vals);
   
%% population plot
M.fig.ax(population_ind) = subplot(nrows,ncols,population_locs);
M.fig.ax(population_ind).Title.String = 'Population Numbers';
M.fig.ax(population_ind).NextPlot = 'add';
M.fig.ax(population_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(population_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(population_ind).XLim = xx;
M.fig.ax(population_ind).Color = [0 0 0];

M.fig.population_plots(1) = plot(M.t,M.setup.N0,'w-','LineWidth',2,'DisplayName','N_T');
M.fig.population_plots(2) = plot(M.t,M.setup.NI0,'color',[0 0.4470 0.7410],'LineWidth',2,'DisplayName','N_I');

legend(M.fig.ax(population_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% subpopulations proportion plot
M.fig.ax(subpop_ind) = subplot(nrows,ncols,subpop_locs);
M.fig.ax(subpop_ind).Title.String = 'Sub-Population Proportions';
M.fig.ax(subpop_ind).NextPlot = 'add';
M.fig.ax(subpop_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(subpop_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(subpop_ind).XLim = xx;
M.fig.ax(subpop_ind).YLim = [0 1];
M.fig.ax(subpop_ind).Color = [0 0 0];

M.fig.subpop_plots(1) = plot(M.t,M.setup.prop_ha0,'w','LineWidth',2,'DisplayName','HA');
M.fig.subpop_plots(2) = plot(M.t,M.setup.prop_mut0,'r','LineWidth',2,'DisplayName','FGFR3 Muts');
M.fig.subpop_plots(3) = plot(M.t,M.setup.NI0/(M.setup.NI0+M.setup.N0),'color',[0 0.4470 0.7410],'LineWidth',2,'DisplayName','CTLs of All Cells');

legend(M.fig.ax(subpop_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% event plot
M.fig.ax(event_ind) = subplot(nrows,ncols,event_locs);
M.fig.ax(event_ind).Title.String = 'Event Rates in Update';
M.fig.ax(event_ind).XLabel.String = 'Days';
M.fig.ax(event_ind).NextPlot = 'add';
M.fig.ax(event_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(event_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(event_ind).XLim = xx;
M.fig.ax(event_ind).YScale = 'log';
M.fig.ax(event_ind).Color = [0 0 0];

M.fig.events_plots(1) = plot(M.t,0,'w','LineWidth',2,'DisplayName','Tumor Apoptosis');
M.fig.events_plots(2) = plot(M.t,0,'w:','LineWidth',2,'DisplayName','Tumor Contact Inhibitions');
M.fig.events_plots(3) = plot(M.t,0,'b','LineWidth',2,'DisplayName','Immune Recruit');
M.fig.events_plots(4) = plot(M.t,0,'g','LineWidth',2,'DisplayName','Immune Prolif');
M.fig.events_plots(5) = plot(M.t,0,'m','LineWidth',2,'DisplayName','Immune Contact Inhibitions');
M.fig.events_plots(6) = plot(M.t,0,'Color',[0.3,0.3,0.9],'LineStyle','-','LineWidth',2,'DisplayName','Immune AICD');

legend(M.fig.ax(event_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% interaction plot
M.fig.ax(interaction_ind) = subplot(nrows,ncols,interaction_locs);
M.fig.ax(interaction_ind).Title.String = 'Interaction Rates in Update';
M.fig.ax(interaction_ind).XLabel.String = 'Days';
M.fig.ax(interaction_ind).NextPlot = 'add';
M.fig.ax(interaction_ind).XLim = xx;
M.fig.ax(interaction_ind).YScale = 'log';
M.fig.ax(interaction_ind).Color = M.fig.ax(interaction_ind).Parent.Color;
M.fig.ax(interaction_ind).Color = [0 0 0];


% M.fig.interaction_plots(1) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','HA Mut IC');
M.fig.interaction_plots(2) = plot(M.t,0,'Color',tumor_colors(3,:),'LineStyle',':','LineWidth',2,'DisplayName','HA Non-Mut IC');
% M.fig.interaction_plots(3) = plot(M.t,0,'Color',tumor_colors(2,:),'LineStyle','-','LineWidth',2,'DisplayName','LA Mut IC');
M.fig.interaction_plots(4) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle',':','LineWidth',2,'DisplayName','LA Non-Mut IC');
M.fig.interaction_plots(5) = plot(M.t,0,'Color',unengaged_color,'LineWidth',2,'DisplayName','Unengaged Deactivations');
M.fig.interaction_plots(6) = plot(M.t,0,'Color',immune_colors(2,:),'LineWidth',2,'DisplayName','Slow-killing Deactivations');
M.fig.interaction_plots(7) = plot(M.t,0,'Color',immune_colors(1,:),'LineWidth',2,'DisplayName','Fast-killing Deactivations');

legend(M.fig.ax(interaction_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% receptor plot
M.fig.ax(receptor_ind) = subplot(nrows,ncols,receptor_locs);
M.fig.ax(receptor_ind).Title.String = 'Receptor proportions';
M.fig.ax(receptor_ind).NextPlot = 'add';
M.fig.ax(receptor_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(receptor_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(receptor_ind).XLim = xx;
M.fig.ax(receptor_ind).Color = [0 0 0];

M.fig.receptor_plots(1) = plot(M.t,0,'w','LineWidth',2,'DisplayName','\phi_D');
M.fig.receptor_plots(2) = plot(M.t,0,'g','LineWidth',2,'DisplayName','Checkpoint');

legend(M.fig.ax(receptor_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% drug plots
M.fig.ax(drug_ind) = subplot(nrows,ncols,drug_locs);
M.fig.ax(drug_ind).Title.String = 'Free Drug Concentrations';
M.fig.ax(drug_ind).NextPlot = 'add';
M.fig.ax(drug_ind).XLim = xx;
M.fig.ax(drug_ind).YScale = 'log';
M.fig.ax(drug_ind).YAxis.TickLabelFormat = '%2.1g';
M.fig.ax(drug_ind).Color = [0 0 0];

if M.fgfr3.n_regions==1
    c = spring(2);
    M.fig.afgfr3_plots(1) = plot(M.t,0,'Color',c(2,:),'LineStyle','-','LineWidth',2,'DisplayName','Free aFGFR3 (R1)');
else
    c = spring(3);
    M.fig.afgfr3_plots(1) = plot(M.t,0,'Color',c(2,:),'LineStyle','-','LineWidth',2,'DisplayName','Free aFGFR3 (R1)');
    M.fig.afgfr3_plots(2) = plot(M.t,0,'Color',c(3,:),'LineStyle','-','LineWidth',2,'DisplayName',sprintf('Free aFGFR3 (R%d)',M.fgfr3.n_regions));
end

M.fig.afgfr3_plots(end+1) = plot(M.t,0,'Color',c(1,:),'LineStyle','-','LineWidth',2,'DisplayName','aFGFR3 in Circ');

if M.checkpoint.n_regions==1
    c = winter(2);
    M.fig.apd1_plots(1) = plot(M.t,0,'Color',c(2,:),'LineStyle','-.','LineWidth',2,'DisplayName','Free aPD1 (R1)');
else
    c = winter(3);
    M.fig.apd1_plots(1) = plot(M.t,0,'Color',c(2,:),'LineStyle','-.','LineWidth',2,'DisplayName','Free aPD1 (R1)');
    M.fig.apd1_plots(2) = plot(M.t,0,'Color',c(3,:),'LineStyle','-.','LineWidth',2,'DisplayName',sprintf('Free aPD1 (R%d)',M.checkpoint.n_regions));
end

M.fig.apd1_plots(end+1) = plot(M.t,0,'Color',c(1,:),'LineStyle','-.','LineWidth',2,'DisplayName','aPD1 in Circ');

legend(M.fig.ax(drug_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)


%% not time series
%% tumor probabilities
M.fig.ax(tum_prob_ind) = subplot(nrows,ncols,tum_probs_locs);
max_tum_prob = sum(1-exp(-M.pars.prolif_rate*M.pars.max_dt));
M.fig.ax(tum_prob_ind).YLim = [0 max_tum_prob];
M.fig.ax(tum_prob_ind).Title.String = 'Tum Outcome Probabilities';
M.fig.ax(tum_prob_ind).NextPlot = 'replacechildren';

M.fig.tum_prob_bar = bar(NaN,NaN(2,1),1,'stacked','EdgeAlpha',0);

events = {'prolif','apop'};
order = [2,1];
legend(M.fig.ax(tum_prob_ind),flip(M.fig.tum_prob_bar),flip(events(order)),'AutoUpdate','off','Location','northwest')

%% immune probabilities
max_imm_prob = min(1,sum(1-exp(-[M.immune_pars.move_rate;M.immune_pars.conjugation_rate;M.immune_pars.apop_rate;M.immune_pars.prolif_rate;0;M.immune_pars.aicd_rate]*M.pars.max_dt_Imm)));
M.fig.ax(imm_prob_ind) = subplot(nrows,ncols,imm_probs_locs);
M.fig.ax(imm_prob_ind).XTick = [];
M.fig.ax(imm_prob_ind).Title.String = 'Imm Outcome Probabilities';
M.fig.ax(imm_prob_ind).YLim = [0 max_imm_prob];
M.fig.ax(imm_prob_ind).NextPlot = 'replacechildren';

M.fig.imm_prob_bar = bar(NaN,NaN(6,1),1,'stacked','EdgeAlpha',0,'Facecolor','flat');

M.fig.imm_events = {'prolif','apop','move','kill','deactivation','AICD'};
M.fig.imm_prob_order = 1:length(M.fig.imm_events);
M.fig.imm_prob_colors = lines(6);
legend(M.fig.ax(imm_prob_ind),flip(M.fig.imm_prob_bar),flip(M.fig.imm_events(M.fig.imm_prob_order)),'AutoUpdate','off','Location','northwest')


drawnow