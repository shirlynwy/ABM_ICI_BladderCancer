function plotFunction_EndStep(M)

%% cell locations
fast_log = M.immunes(:,M.I.type)==2;
slow_log = M.immunes(:,M.I.type)==1;
unengaged_log = M.immunes(:,M.I.type)==0;
deactivated_log = M.immunes(:,M.I.type)==-1;
if M.plot_pars.plot_location

    data_names = ["XData","YData","ZData"];

    sc_data = cell(8,1);
    sc_data{1} = M.tumors(~M.ha_log & ~M.mut_log,M.I.subs)-M.grid.center;
%     sc_data{2} = M.tumors(~M.ha_log & M.mut_log,M.I.subs)-M.grid.center;
    sc_data{3} = M.tumors(M.ha_log & ~M.mut_log,M.I.subs)-M.grid.center;
%     sc_data{4} = M.tumors(M.ha_log & M.mut_log,M.I.subs)-M.grid.center;
    sc_data{5} = M.immunes(fast_log,M.I.subs)-M.grid.center;
    sc_data{6} = M.immunes(slow_log,M.I.subs)-M.grid.center;
    sc_data{7} = M.immunes(unengaged_log,M.I.subs)-M.grid.center;
    sc_data{8} = M.immunes(deactivated_log,M.I.subs)-M.grid.center;

    if M.plot_pars.use_carveout
        th = (M.t*0.2 - 0.05)*pi; % shift the carve out so it is not head on
        for axi =[1,3, 5:8]%1:8
            I = ( sc_data{axi}(:,3) <= 0 ) | ( sqrt(2)*(sc_data{axi}(:,1:2)*[cos(th);sin(th)]) < hypot(sc_data{axi}(:,1),sc_data{axi}(:,2)) );
            sc_data{axi} = sc_data{axi}(I,:);
        end
    end
    for axi = [1,3, 5:8]%1:8
        for di = 1:3
            M.fig.scatter_plots(axi).(data_names(di)) = sc_data{axi}(:,di);
        end
    end

%     M.fig.scatter_plots(1).XData = M.tumors(M.ha_log,M.I.subs(1))-M.grid.center(1);
%     M.fig.scatter_plots(1).YData = M.tumors(M.ha_log,M.I.subs(2))-M.grid.center(2);
%     M.fig.scatter_plots(1).ZData = M.tumors(M.ha_log,M.I.subs(3))-M.grid.center(3);
%     
%     M.fig.scatter_plots(2).XData = M.tumors(~M.ha_log,M.I.subs(1))-M.grid.center(1);
%     M.fig.scatter_plots(2).YData = M.tumors(~M.ha_log,M.I.subs(2))-M.grid.center(2);
%     M.fig.scatter_plots(2).ZData = M.tumors(~M.ha_log,M.I.subs(3))-M.grid.center(3);
%   
%     M.fig.scatter_plots(3).XData = M.immunes(fast_log,M.I.subs(1))-M.grid.center(1);
%     M.fig.scatter_plots(3).YData = M.immunes(fast_log,M.I.subs(2))-M.grid.center(2);
%     M.fig.scatter_plots(3).ZData = M.immunes(fast_log,M.I.subs(3))-M.grid.center(3);
%     
%     M.fig.scatter_plots(4).XData = M.immunes(slow_log,M.I.subs(1))-M.grid.center(1);
%     M.fig.scatter_plots(4).YData = M.immunes(slow_log,M.I.subs(2))-M.grid.center(2);
%     M.fig.scatter_plots(4).ZData = M.immunes(slow_log,M.I.subs(3))-M.grid.center(3);
%     
%     M.fig.scatter_plots(5).XData = M.immunes(unengaged_log,M.I.subs(1))-M.grid.center(1);
%     M.fig.scatter_plots(5).YData = M.immunes(unengaged_log,M.I.subs(2))-M.grid.center(2);
%     M.fig.scatter_plots(5).ZData = M.immunes(unengaged_log,M.I.subs(3))-M.grid.center(3);
%     
%     M.fig.scatter_plots(6).XData = M.immunes(deactivated_log,M.I.subs(1))-M.grid.center(1);
%     M.fig.scatter_plots(6).YData = M.immunes(deactivated_log,M.I.subs(2))-M.grid.center(2);
%     M.fig.scatter_plots(6).ZData = M.immunes(deactivated_log,M.I.subs(3))-M.grid.center(3);

    view(M.fig.ax(M.fig.scatter_ind), M.plot_pars.m*[cos(M.t*0.2*pi)*0.8660,sin(M.t*0.2*pi)*0.8660,0.5]) % full rotation every 10 days
%     view(M.fig.ax(M.fig.scatter_ind),[36*M.t 30 40])
end

%% slice
z_mid_log_tumor = M.tumors(:,M.I.subs(3))==M.grid.center(3);
z_mid_log_immune = M.immunes(:,M.I.subs(3))==M.grid.center(3);
M.fig.cell_slice_plot(1).XData = M.tumors(~M.ha_log & ~M.mut_log & z_mid_log_tumor,M.I.subs(1));
M.fig.cell_slice_plot(1).YData = M.tumors(~M.ha_log & ~M.mut_log & z_mid_log_tumor,M.I.subs(2));
M.fig.cell_slice_plot(1).ZData = ones(length(M.fig.cell_slice_plot(1).XData),1);

% M.fig.cell_slice_plot(2).XData = M.tumors(~M.ha_log & M.mut_log & z_mid_log_tumor,M.I.subs(1));
% M.fig.cell_slice_plot(2).YData = M.tumors(~M.ha_log & M.mut_log & z_mid_log_tumor,M.I.subs(2));
% M.fig.cell_slice_plot(2).ZData = ones(length(M.fig.cell_slice_plot(2).XData),1);

M.fig.cell_slice_plot(3).XData = M.tumors(M.ha_log & ~M.mut_log & z_mid_log_tumor,M.I.subs(1));
M.fig.cell_slice_plot(3).YData = M.tumors(M.ha_log & ~M.mut_log & z_mid_log_tumor,M.I.subs(2));
M.fig.cell_slice_plot(3).ZData = ones(length(M.fig.cell_slice_plot(3).XData),1);

% M.fig.cell_slice_plot(4).XData = M.tumors(M.ha_log & M.mut_log & z_mid_log_tumor,M.I.subs(1));
% M.fig.cell_slice_plot(4).YData = M.tumors(M.ha_log & M.mut_log & z_mid_log_tumor,M.I.subs(2));
% M.fig.cell_slice_plot(4).ZData = ones(length(M.fig.cell_slice_plot(4).XData),1);

M.fig.cell_slice_plot(5).XData = M.immunes(fast_log & z_mid_log_immune,M.I.subs(1));
M.fig.cell_slice_plot(5).YData = M.immunes(fast_log & z_mid_log_immune,M.I.subs(2));
M.fig.cell_slice_plot(5).ZData = ones(length(M.fig.cell_slice_plot(5).XData),1);

M.fig.cell_slice_plot(6).XData = M.immunes(slow_log & z_mid_log_immune,M.I.subs(1));
M.fig.cell_slice_plot(6).YData = M.immunes(slow_log & z_mid_log_immune,M.I.subs(2));
M.fig.cell_slice_plot(6).ZData = ones(length(M.fig.cell_slice_plot(6).XData),1);

M.fig.cell_slice_plot(7).XData = M.immunes(unengaged_log & z_mid_log_immune,M.I.subs(1));
M.fig.cell_slice_plot(7).YData = M.immunes(unengaged_log & z_mid_log_immune,M.I.subs(2));
M.fig.cell_slice_plot(7).ZData = ones(length(M.fig.cell_slice_plot(7).XData),1);

M.fig.cell_slice_plot(8).XData = M.immunes(deactivated_log & z_mid_log_immune,M.I.subs(1));
M.fig.cell_slice_plot(8).YData = M.immunes(deactivated_log & z_mid_log_immune,M.I.subs(2));
M.fig.cell_slice_plot(8).ZData = ones(length(M.fig.cell_slice_plot(8).XData),1);

M.fig.cell_slice_plot(9).ZData = M.immune_stimulatory_factor.concentration(:,:,M.grid.center(3))';

view(M.fig.ax(M.fig.cell_slice_ind),2);

%% projections
if M.NT>1
    [N,~,~] = histcounts2(M.tumors(:,M.I.subs(1)),M.tumors(:,M.I.subs(2)),(M.grid.xx(1)-.5):(M.grid.xx(end)+.5),(M.grid.yy(1)-.5):(M.grid.yy(end)+.5),'Normalization','countdensity');
    if all(size(N)>1)
        M.fig.tum_density_plot.ZData = N';
    else
        M.fig.tum_density_plot.ZData = [];        
    end
end


if M.NI>1 && M.NT>1
    [N,~,~] = histcounts2(M.immunes(~deactivated_log,M.I.subs(1)),M.immunes(~deactivated_log,M.I.subs(2)),(M.grid.xx(1)-.5):(M.grid.xx(end)+.5),(M.grid.yy(1)-.5):(M.grid.yy(end)+.5),'Normalization','countdensity');
    if all(size(N)>1)
        M.fig.imm_density_plot.ZData = N';
    else
        M.fig.imm_density_plot.ZData = [];        
    end
end

M.fig.isf_density_plot.ZData = mean(M.immune_stimulatory_factor.concentration,3)';
M.fig.afgfr3_density_plot.ZData = M.fgfr3.concentration(M.fgfr3.regions(:,:,M.grid.center(3)))';
M.fig.apd1_density_plot.ZData = M.checkpoint.aPD1(M.checkpoint.regions(:,:,M.grid.center(3)))';

%% tracked quantities
%% population plot
M.fig.population_plots(1).XData(end+1) = M.t;
M.fig.population_plots(1).YData(end+1) = M.NT;
M.fig.population_plots(2).XData(end+1) = M.t;
M.fig.population_plots(2).YData(end+1) = M.NI;

%% subpopulations proportion plot
M.fig.subpop_plots(1).XData(end+1) = M.t;
M.fig.subpop_plots(1).YData(end+1) = sum(M.tumors(:,M.I.type)==1)./M.NT; % HA
M.fig.subpop_plots(2).XData(end+1) = M.t;
M.fig.subpop_plots(2).YData(end+1) = sum(M.tumors(:,M.I.tumor_mut))./M.NT; % mut
M.fig.subpop_plots(3).XData(end+1) = M.t;
M.fig.subpop_plots(3).YData(end+1) = M.NI/(M.NI+M.NT); % immune

%% receptor plot
M.fig.receptor_plots(1).XData(end+1) = M.t;
M.fig.receptor_plots(1).YData(end+1) = M.tracked.phiD_mean(M.i) / M.fgfr3.gammaT;

%% drug plot
M.fig.afgfr3_plots(1).XData(end+1) = M.t;
if M.fgfr3.concentration(1) < 1e-2
    M.fig.afgfr3_plots(1).YData(end+1) = 0;
else
    M.fig.afgfr3_plots(1).YData(end+1) = M.fgfr3.concentration(1);
end
if M.fgfr3.n_regions>1
    M.fig.afgfr3_plots(2).XData(end+1) = M.t;
    if M.fgfr3.concentration(end) < 1e-2
        M.fig.afgfr3_plots(2).YData(end+1) = 0;
    else
        M.fig.afgfr3_plots(2).YData(end+1) = M.fgfr3.concentration(end);
    end
end
M.fig.afgfr3_plots(end).XData(end+1) = M.t;
if M.fgfr3.circ < 1e-2
    M.fig.afgfr3_plots(end).YData(end+1) = 0;
else
    M.fig.afgfr3_plots(end).YData(end+1) = M.fgfr3.circ;
end

M.fig.apd1_plots(1).XData(end+1) = M.t;
if M.checkpoint.aPD1(1) < 1e-2
    M.fig.apd1_plots(1).YData(end+1) = 0;
else
    M.fig.apd1_plots(1).YData(end+1) = M.checkpoint.aPD1(1);
end
if M.fgfr3.n_regions>1
    M.fig.apd1_plots(2).XData(end+1) = M.t;
    if M.checkpoint.aPD1(end) < 1e-2
        M.fig.apd1_plots(2).YData(end+1) = 0;
    else
        M.fig.apd1_plots(2).YData(end+1) = M.checkpoint.aPD1(end);
    end
end
M.fig.apd1_plots(end).XData(end+1) = M.t;
if M.checkpoint.aPD1_circulation < 1e-2
    M.fig.apd1_plots(end).YData(end+1) = 0;
else
    M.fig.apd1_plots(end).YData(end+1) = M.checkpoint.aPD1_circulation;
end

%% event plots
M.fig.events_plots(1).XData(end+1) = M.t;
M.fig.events_plots(1).YData(end+1) = sum(M.tracked.tum_apop(M.i,:,:),2:3) / (M.NT*M.dt);
M.fig.events_plots(2).XData(end+1) = M.t;
M.fig.events_plots(2).YData(end+1) = M.tracked.tum_contact_inhibition(M.i) / (M.NT*M.dt);
M.fig.events_plots(3).XData(end+1) = M.t;
M.fig.events_plots(3).YData(end+1) = M.tracked.imm_recruit(M.i) / M.dt;
M.fig.events_plots(4).XData(end+1) = M.t;
M.fig.events_plots(4).YData(end+1) = M.tracked.imm_prolif(M.i) / (M.NI*M.dt);
M.fig.events_plots(5).XData(end+1) = M.t;
M.fig.events_plots(5).YData(end+1) = M.tracked.imm_contact_inhibition(M.i) / (M.NI*M.dt);
M.fig.events_plots(6).XData(end+1) = M.t;
M.fig.events_plots(6).YData(end+1) = M.tracked.imm_aicd(M.i) / (M.NI*M.dt);

%% interaction plots
% M.fig.interaction_plots(1).XData(end+1) = M.t;
% M.fig.interaction_plots(1).YData(end+1) = M.tracked.imm_cleared(M.i,2,2) / (sum(M.ha_log & M.mut_log)*M.dt);
M.fig.interaction_plots(2).XData(end+1) = M.t;
M.fig.interaction_plots(2).YData(end+1) = M.tracked.imm_cleared(M.i,1,2) / (sum(M.ha_log & ~M.mut_log)*M.dt);
% M.fig.interaction_plots(3).XData(end+1) = M.t;
% M.fig.interaction_plots(3).YData(end+1) = M.tracked.imm_cleared(M.i,2,1) / (sum(~M.ha_log & M.mut_log)*M.dt);
M.fig.interaction_plots(4).XData(end+1) = M.t;
M.fig.interaction_plots(4).YData(end+1) = M.tracked.imm_cleared(M.i,1,1) / (sum(~M.ha_log & ~M.mut_log)*M.dt);
M.fig.interaction_plots(5).XData(end+1) = M.t;
M.fig.interaction_plots(5).YData(end+1) = M.tracked.deactivations(M.i,1) / (sum(unengaged_log)*M.dt);
M.fig.interaction_plots(6).XData(end+1) = M.t;
M.fig.interaction_plots(6).YData(end+1) = M.tracked.deactivations(M.i,2) / (sum(slow_log)*M.dt);
M.fig.interaction_plots(7).XData(end+1) = M.t;
M.fig.interaction_plots(7).YData(end+1) = M.tracked.deactivations(M.i,3) / (sum(fast_log)*M.dt);

drawnow